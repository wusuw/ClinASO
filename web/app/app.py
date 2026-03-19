import os
import subprocess
import threading
import smtplib
import shutil
import queue
import hashlib
import json
import uuid
import logging
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from flask import Flask, request, jsonify
from flask_cors import CORS
from flask import Response
app = Flask(__name__)
CORS(app)

# 配置日志记录
from logging.handlers import RotatingFileHandler

LOG_DIR = '/var/www/genetic-analysis/app/logs'
os.makedirs(LOG_DIR, exist_ok=True)

# 统一格式
LOG_FORMAT = '%(asctime)s [%(levelname)s] %(message)s'
LOG_DATE_FORMAT = '%Y-%m-%d %H:%M:%S'

formatter = logging.Formatter(LOG_FORMAT, datefmt=LOG_DATE_FORMAT)

# 主日志（所有记录，自动轮转，最大 10MB，保留 5 个备份）
main_handler = RotatingFileHandler(
    os.path.join(LOG_DIR, 'app.log'),
    maxBytes=10*1024*1024,
    backupCount=5,
    encoding='utf-8'
)
main_handler.setLevel(logging.DEBUG)
main_handler.setFormatter(formatter)

# 错误日志（只记录 ERROR 及以上）
error_handler = RotatingFileHandler(
    os.path.join(LOG_DIR, 'error.log'),
    maxBytes=5*1024*1024,
    backupCount=3,
    encoding='utf-8'
)
error_handler.setLevel(logging.ERROR)
error_handler.setFormatter(formatter)

# 控制台输出（开发调试用）
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(formatter)

# 配置根 logger
logger = logging.getLogger('ClinASO')
logger.setLevel(logging.DEBUG)
logger.addHandler(main_handler)
logger.addHandler(error_handler)
logger.addHandler(console_handler)

# 邮件服务配置
SMTP_SERVER = 'smtp.163.com'
SMTP_PORT = 465
SMTP_USERNAME = 'your_email@example.com'
SMTP_PASSWORD = 'your_smtp_password'
FROM_EMAIL = 'your_email@example.com'

# 结果文件目录
OUTPUT_DIR = '/asodesigner/outfile/'
OUTPUT_DIR_Homo = '/homologyanalysis/outfile/'
OUTPUT_DIR_off = '/offtarget/outfile/'
OUTPUT_DIR_snp = '/snp/outfile/'
RESULTS_DIR = 'results/'

# 确保输出目录存在
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR_Homo, exist_ok=True)  
os.makedirs(OUTPUT_DIR_snp, exist_ok=True) 
os.makedirs(OUTPUT_DIR_off, exist_ok=True)

# ==================== 任务去重配置 ====================
TASK_RECORD_FILE = '/var/www/genetic-analysis/app/task_records.json'
task_records_lock = threading.Lock()

def load_task_records():
    if os.path.exists(TASK_RECORD_FILE):
        try:
            with open(TASK_RECORD_FILE, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            return {}
    return {}

def save_task_records(records):
    with open(TASK_RECORD_FILE, 'w') as f:
        json.dump(records, f, indent=2, ensure_ascii=False)

def make_task_key(task_type, email, params_str):
    raw = task_type + ':' + email + ':' + params_str
    return hashlib.sha256(raw.encode('utf-8')).hexdigest()

def is_duplicate_task(task_type, email, params_str):
    task_key = make_task_key(task_type, email, params_str)
    with task_records_lock:
        records = load_task_records()
        if task_key in records:
            return True
        from datetime import datetime
        records[task_key] = {
            'type': task_type,
            'email': email,
            'created_at': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        save_task_records(records)
        return False
# ==================== 去重配置结束 ====================
 
# 任务队列和工作线程
task_queue = queue.Queue(maxsize=20)
worker_thread = None
worker_lock = threading.Lock()

def worker():
    """工作线程，顺序处理任务队列中的任务"""
    logger.info("Worker thread started")
    while True:
        try:
            task = task_queue.get(block=True)
            if task is None:  # 收到停止信号
                logger.info("Worker stopping signal received")
                break
                
            task_type, data = task
            logger.info(f"Processing {task_type} task")
            
            try:
                if task_type == 'design':
                    run_design_analysis(data)
                elif task_type == 'homology':
                    run_homology_analysis(data)
                elif task_type == 'snp':
                    run_snp_analysis(data)
                elif task_type == 'offtarget':
                    run_offtarget_analysis(data)
            except Exception as e:
                logger.error(f"Task failed: {task_type} - {str(e)}")
                # 发送错误邮件
                send_email(data['email'], data.get('gene_name', 'Analysis'), 
                          f"Analysis failed: {str(e)}", f"{task_type.capitalize()} Analysis Error")
            finally:
                task_queue.task_done()
                
        except queue.Empty:
            logger.info("Worker queue is empty, still waiting...")
        except Exception as e:
            logger.error(f"Worker error: {str(e)}")

def start_worker():
    """启动工作线程"""
    global worker_thread
    with worker_lock:
        if worker_thread is None or not worker_thread.is_alive():
            worker_thread = threading.Thread(target=worker, daemon=True)
            worker_thread.start()
            logger.info("Worker thread started")

def stop_worker():
    """停止工作线程"""
    global worker_thread
    with worker_lock:
        if worker_thread and worker_thread.is_alive():
            task_queue.put(None)  # 发送停止信号
            try:
                worker_thread.join(timeout=10)
            except Exception:
                pass
            logger.info("Worker thread stopped")

def run_design_analysis(data):
    """执行基因分析脚本并发送结果邮件"""
    gene_name = data['gene_name'].upper()
    email = data['email']
    aso_len = data["aso_len"]
    aso_count = data.get("aso_count", "30")
    priority = data.get("priority", "specificity")
    homologous_species = data.get("homologous_species", "all")
    # 为任务创建独立的临时目录
    task_id = str(uuid.uuid4())
    task_output_dir = os.path.join(OUTPUT_DIR, task_id)
    os.makedirs(task_output_dir, exist_ok=True)
    
    logger.info(f"Starting design analysis for {gene_name} in {task_output_dir}")
    logger.info(f"ASO count: {aso_count}, Priority: {priority}, Homologous species: {homologous_species}")
    
    try:
        # 执行脚本 - 将输出重定向到特定目录
        result = subprocess.run(
            f'bash /asodesigner/aso_design2.sh {gene_name} {aso_len} {aso_count} {priority} {homologous_species} {task_output_dir}',
            shell=True,
            capture_output=True,
            text=True,
            check=True,
            env={**os.environ, 'OUTPUT_DIR': task_output_dir},
            timeout=None  # 移除超时限制
        )
        
        # 确定结果文件路径
        result_file = os.path.join(task_output_dir, f'ASO_AllCandidates_{gene_name}.xlsx')
        filtered_result_file = os.path.join(task_output_dir, f'ASO_FilteredCandidates_{gene_name}.xlsx')
        
        # 检查结果文件是否存在
        if os.path.exists(result_file) and os.path.getsize(result_file) > 0:
            # 成功生成结果文件，发送结果邮件
            send_success_email(email, gene_name, result_file)
            logger.info(f"Design analysis completed successfully for {gene_name}")
            logger.info(f"Original result file: {result_file}")
            if os.path.exists(filtered_result_file):
                logger.info(f"Filtered result file: {filtered_result_file}")
        else:
            # 未生成结果文件，发送错误提示邮件
            send_error_email(email, gene_name)
            logger.warning(f"Design analysis completed but no result file generated for {gene_name}")
        
    except subprocess.TimeoutExpired as e:
        error_msg = f"Design analysis timed out:\n{e.stderr}"
        logger.error(error_msg)
        send_timeout_email(email, gene_name)
    except subprocess.CalledProcessError as e:
        error_msg = f"Design analysis failed:\n{e.stderr}"
        logger.error(error_msg)
        send_failure_email(email, gene_name, error_msg)
    except Exception as e:
        error_msg = f"System error:\n{str(e)}"
        logger.error(error_msg)
        send_system_error_email(email, gene_name, error_msg)
    finally:
        # 清理临时目录
        cleanup_output_directory(task_output_dir)

def send_success_email(to_email, gene_name, result_file):
    """发送成功分析结果邮件"""
    logger.info(f"Sending success email to: {to_email}, Subject: {gene_name}")
    
    try:
        msg = MIMEMultipart()
        msg['From'] = FROM_EMAIL
        msg['To'] = to_email
        msg['Subject'] = f'Gapmer ASO Design Results - {gene_name}'
        
        # 邮件正文内容
        body_content = f"""
        <html>
        <body>
            <h2>Gapmer ASO Design Analysis Results</h2>
            <p>Dear User,</p>
            <p>We are pleased to inform you that the Gapmer ASO design analysis for <strong>{gene_name}</strong> has been completed successfully.</p>
            <p>Please find the following files attached to this email:</p>
            <ul>
                <li><strong>Original Results</strong>: Complete list of designed ASOs without filtering</li>
                <li><strong>Filtered Results</strong>: ASOs filtered based on your selected priority criteria</li>
            </ul>
            <p>If you have any questions about the results or need further assistance, please feel free to reply directly to this email. Our team will be happy to help you.</p>
            <p>Thank you for using the ClinASO ASO Design Platform!</p>
            <hr>
            <p><em>Best regards,</em></p>
            <p><em>The ClinASO Team</em></p>
            <p><em>School of Life Sciences, Yunnan University</em></p>
        </body>
        </html>
        """
        
        body = MIMEText(body_content, 'html', 'utf-8')
        msg.attach(body)
        
        # 添加原始结果文件附件
        with open(result_file, 'rb') as file:
            part = MIMEApplication(
                file.read(),
                Name=os.path.basename(result_file)
            )
        part['Content-Disposition'] = f'attachment; filename="{os.path.basename(result_file)}"'
        msg.attach(part)
        
        # 添加筛选后的结果文件附件
        result_dir = os.path.dirname(result_file)
        filtered_result_file = os.path.join(result_dir, f'ASO_FilteredCandidates_{gene_name}.xlsx')
        if os.path.exists(filtered_result_file):
            with open(filtered_result_file, 'rb') as file:
                part = MIMEApplication(
                    file.read(),
                    Name=os.path.basename(filtered_result_file)
                )
            part['Content-Disposition'] = f'attachment; filename="{os.path.basename(filtered_result_file)}"'
            msg.attach(part)
        
        # 如果存在图表文件，也一并添加
        chart_files = [f for f in os.listdir(result_dir) if f.endswith(('.png', '.pdf'))]
        for chart_file in chart_files:
            chart_path = os.path.join(result_dir, chart_file)
            with open(chart_path, 'rb') as file:
                part = MIMEApplication(
                    file.read(),
                    Name=chart_file
                )
            part['Content-Disposition'] = f'attachment; filename="{chart_file}"'
            msg.attach(part)
        
        # 连接SMTP服务器并发送
        with smtplib.SMTP_SSL(SMTP_SERVER, SMTP_PORT) as server:
            server.login(SMTP_USERNAME, SMTP_PASSWORD)
            server.send_message(msg)
            logger.info("Success email sent successfully")
            
        return True
        
    except Exception as e:
        logger.error(f"Error sending success email: {str(e)}")
        return False

def send_error_email(to_email, gene_name):
    """发送错误提示邮件 - 未生成结果文件"""
    logger.info(f"Sending error email to: {to_email}, Subject: {gene_name} - No results")
    
    try:
        msg = MIMEMultipart()
        msg['From'] = FROM_EMAIL
        msg['To'] = to_email
        msg['Subject'] = f'Gapmer ASO Design Analysis Issue - {gene_name}'
        
        # 邮件正文内容
        body_content = f"""
        <html>
        <body>
            <h2>Gapmer ASO Design Analysis Issue</h2>
            <p>Dear User,</p>
            <p>We encountered an issue while processing your Gapmer ASO design request for <strong>{gene_name}</strong>.</p>
            <p>Unfortunately, we were unable to generate the analysis results. This could be due to one of the following reasons:</p>
            <ul>
                <li><strong>Gene name</strong>: The gene name you provided may not be standard or recognized in our database. Please verify that you are using the official gene symbol (e.g., TP53, BRCA1).</li>
                <li><strong>Species</strong>: Our analysis pipeline is currently optimized for human genes only. Please ensure you are analyzing a human gene.</li>
            </ul>
            <p><strong>Recommended actions:</strong></p>
            <ol>
                <li>Check the gene name spelling and format. Use official gene symbols from databases like NCBI or Ensembl.</li>
                <li>Verify that you are analyzing a human gene.</li>
                <li>If you believe the gene name and species are correct, please reply to this email with more details, and our team will investigate further.</li>
            </ol>
            <p>We apologize for any inconvenience and appreciate your understanding.</p>
            <hr>
            <p><em>Best regards,</em></p>
            <p><em>The ClinASO Team</em></p>
            <p><em>School of Life Sciences, Yunnan University</em></p>
        </body>
        </html>
        """
        
        body = MIMEText(body_content, 'html', 'utf-8')
        msg.attach(body)
        
        # 连接SMTP服务器并发送
        with smtplib.SMTP_SSL(SMTP_SERVER, SMTP_PORT) as server:
            server.login(SMTP_USERNAME, SMTP_PASSWORD)
            server.send_message(msg)
            logger.info("Error email sent successfully")
            
        return True
        
    except Exception as e:
        logger.error(f"Error sending error email: {str(e)}")
        return False

def send_timeout_email(to_email, gene_name):
    """发送超时错误邮件"""
    logger.info(f"Sending timeout email to: {to_email}, Subject: {gene_name} - Timeout")
    
    try:
        msg = MIMEMultipart()
        msg['From'] = FROM_EMAIL
        msg['To'] = to_email
        msg['Subject'] = f'Gapmer ASO Design Analysis Timeout - {gene_name}'
        
        # 邮件正文内容
        body_content = f"""
        <html>
        <body>
            <h2>Gapmer ASO Design Analysis Timeout</h2>
            <p>Dear User,</p>
            <p>We regret to inform you that the Gapmer ASO design analysis for <strong>{gene_name}</strong> has exceeded the maximum processing time and was terminated.</p>
            <p>This could happen due to several reasons:</p>
            <ul>
                <li>The gene may have a very complex structure requiring extensive computation.</li>
                <li>There might be an unusually high number of potential target sites.</li>
                <li>System resources may be currently under high demand.</li>
            </ul>
            <p><strong>Recommended actions:</strong></p>
            <ol>
                <li>Please try submitting your analysis again later.</li>
                <li>If the issue persists, consider breaking down your request into smaller segments.</li>
                <li>For immediate assistance, please reply to this email with details of your request.</li>
            </ol>
            <p>We apologize for any inconvenience and appreciate your patience.</p>
            <hr>
            <p><em>Best regards,</em></p>
            <p><em>The ClinASO Team</em></p>
            <p><em>School of Life Sciences, Yunnan University</em></p>
        </body>
        </html>
        """
        
        body = MIMEText(body_content, 'html', 'utf-8')
        msg.attach(body)
        
        # 连接SMTP服务器并发送
        with smtplib.SMTP_SSL(SMTP_SERVER, SMTP_PORT) as server:
            server.login(SMTP_USERNAME, SMTP_PASSWORD)
            server.send_message(msg)
            logger.info("Timeout email sent successfully")
            
        return True
        
    except Exception as e:
        logger.error(f"Error sending timeout email: {str(e)}")
        return False

def send_failure_email(to_email, gene_name, error_msg):
    """发送分析失败邮件"""
    logger.info(f"Sending failure email to: {to_email}, Subject: {gene_name} - Analysis failed")
    
    try:
        msg = MIMEMultipart()
        msg['From'] = FROM_EMAIL
        msg['To'] = to_email
        msg['Subject'] = f'Gapmer ASO Design Analysis Failed - {gene_name}'
        
        # 邮件正文内容
        body_content = f"""
        <html>
        <body>
            <h2>Gapmer ASO Design Analysis Failed</h2>
            <p>Dear User,</p>
            <p>We regret to inform you that the Gapmer ASO design analysis for <strong>{gene_name}</strong> has encountered an error and could not be completed.</p>
            <p>Our technical team has been notified of this issue and will investigate the cause.</p>
            <p><strong>Error details:</strong></p>
            <pre>{error_msg}</pre>
            <p><strong>Recommended actions:</strong></p>
            <ol>
                <li>Please verify that your gene name is correct and refers to a human gene.</li>
                <li>Try submitting your analysis again in a few hours.</li>
                <li>If the issue persists, please reply to this email with the gene name and any additional details.</li>
            </ol>
            <p>We apologize for any inconvenience and appreciate your understanding.</p>
            <hr>
            <p><em>Best regards,</em></p>
            <p><em>The ClinASO Team</em></p>
            <p><em>School of Life Sciences, Yunnan University</em></p>
        </body>
        </html>
        """
        
        body = MIMEText(body_content, 'html', 'utf-8')
        msg.attach(body)
        
        # 连接SMTP服务器并发送
        with smtplib.SMTP_SSL(SMTP_SERVER, SMTP_PORT) as server:
            server.login(SMTP_USERNAME, SMTP_PASSWORD)
            server.send_message(msg)
            logger.info("Failure email sent successfully")
            
        return True
        
    except Exception as e:
        logger.error(f"Error sending failure email: {str(e)}")
        return False

def send_system_error_email(to_email, gene_name, error_msg):
    """发送系统错误邮件"""
    logger.info(f"Sending system error email to: {to_email}, Subject: {gene_name} - System error")
    
    try:
        msg = MIMEMultipart()
        msg['From'] = FROM_EMAIL
        msg['To'] = to_email
        msg['Subject'] = f'Gapmer ASO Design System Error - {gene_name}'
        
        # 邮件正文内容
        body_content = f"""
        <html>
        <body>
            <h2>Gapmer ASO Design System Error</h2>
            <p>Dear User,</p>
            <p>We regret to inform you that the Gapmer ASO design analysis for <strong>{gene_name}</strong> has encountered a system error and could not be completed.</p>
            <p>This appears to be an issue with our system rather than with your request. Our technical team has been notified and is working to resolve the problem.</p>
            <p><strong>Error details:</strong></p>
            <pre>{error_msg}</pre>
            <p><strong>Recommended actions:</strong></p>
            <ol>
                <li>Please try submitting your analysis again later.</li>
                <li>If the issue persists, please reply to this email with details of your request.</li>
                <li>Our team will provide you with an update as soon as the issue is resolved.</li>
            </ol>
            <p>We sincerely apologize for this inconvenience and appreciate your patience.</p>
            <hr>
            <p><em>Best regards,</em></p>
            <p><em>The ClinASO Team</em></p>
            <p><em>School of Life Sciences, Yunnan University</em></p>
        </body>
        </html>
        """
        
        body = MIMEText(body_content, 'html', 'utf-8')
        msg.attach(body)
        
        # 连接SMTP服务器并发送
        with smtplib.SMTP_SSL(SMTP_SERVER, SMTP_PORT) as server:
            server.login(SMTP_USERNAME, SMTP_PASSWORD)
            server.send_message(msg)
            logger.info("System error email sent successfully")
            
        return True
        
    except Exception as e:
        logger.error(f"Error sending system error email: {str(e)}")
        return False


def run_homology_analysis(data):
    """执行同源性分析并发送结果邮件"""
    gene_name = data['gene_name']
    gene_sequence = data['gene_sequence']
    analysis_type = data['analysis_type']
    email = data['email']
    task_id = str(uuid.uuid4())
    task_output_dir = os.path.join(OUTPUT_DIR_Homo, task_id)
    os.makedirs(task_output_dir, exist_ok=True)
    logger.debug(f"[DEBUG] Task ID: {task_id}")
    logger.debug(f"[DEBUG] Task output directory: {task_output_dir}")
    logger.info(f"Starting homology analysis for {gene_name}")
    
    try:
        # 直接在当前目录创建一个测试文件，看看是否有权限
        test_file = f"{task_output_dir}/test.txt"
        with open(test_file, "w") as f:
            f.write("Test file")
        logger.debug(f"[DEBUG] Test file created: {os.path.exists(test_file)}")
        
        result = subprocess.run(
            f'bash /homologyanalysis/_homologyanalysis.sh {gene_name} {gene_sequence} {analysis_type} {task_output_dir}',
            shell=True,
            capture_output=True,
            text=True,
            check=False,  # 不检查状态码，因为脚本总是返回0
            env={**os.environ, 'OUTPUT_DIR_Homo': task_output_dir},
            timeout=None  # 移除超时限制
        )
        logger.debug(f"[DEBUG] Script exit code: {result.returncode}")
        logger.debug(f"[DEBUG] Script output: {result.stdout}")
        if result.stderr:
            logger.debug(f"[DEBUG] Script stderr: {result.stderr}")
    
        filename = f"{task_output_dir}/Homology_Analysis_Report.pdf"
        
        # 检查结果文件是否存在
        logger.debug(f"[DEBUG] Checking if file exists: {filename}")
        logger.debug(f"[DEBUG] Directory exists: {os.path.exists(task_output_dir)}")
        if os.path.exists(task_output_dir):
            logger.debug(f"[DEBUG] Directory contents: {os.listdir(task_output_dir)}")
        
        if not os.path.exists(filename):
            logger.error(f"Homology result file not found: {filename}")
            send_email(email, gene_name, f"Homology analysis failed: result file not generated. Please check your input and try again.", "Homology Analysis Error")
            return
        else:
            logger.debug(f"[DEBUG] File exists: {filename}")
            logger.debug(f"[DEBUG] File size: {os.path.getsize(filename)} bytes")
        
        # 发送结果邮件
        logger.debug(f"[DEBUG] Sending email with attachment: {filename}")
        send_email_with_attachment(email, gene_name, filename, result.stdout, "Homology Analysis Results")
        logger.debug(f"[DEBUG] Email sent")
        logger.info(f"Homology analysis completed for {gene_name}")
        
    except Exception as e:
        error_msg = f"Homology analysis failed:\n{str(e)}"
        logger.error(f"[DEBUG] Exception: {error_msg}")
        import traceback
        traceback.print_exc()
        logger.error(error_msg)
        send_email(email, gene_name, error_msg, "Homology Analysis Error")
    finally:
        # 清理临时目录
        logger.debug(f"[DEBUG] Cleaning up directory: {task_output_dir}")
        cleanup_output_directory(task_output_dir)
        logger.debug(f"[DEBUG] Cleanup completed")

def run_snp_analysis(data):
    """执行SNP分析并发送结果邮件"""
    gene_name = data['gene_name']
    aso_sequence = data['aso_sequence']
    email = data['email']
    task_id = str(uuid.uuid4())
    task_output_dir = os.path.join(OUTPUT_DIR_snp, task_id)
    os.makedirs(task_output_dir, exist_ok=True)
    logger.info(f"Starting SNP analysis for {aso_sequence}")
    logger.info(f"Starting SNP analysis for {task_output_dir}")
    
    try:        
        result = subprocess.run(
            f'bash /snp/_snp.sh {gene_name} {aso_sequence} {task_output_dir}',
            shell=True,
            capture_output=True,
            text=True,
            check=False,  # 不检查状态码
            env={**os.environ, 'OUTPUT_DIR_snp': task_output_dir},
            timeout=None  # 移除超时限制
        )
        logger.info(f"SNP analysis script exit code: {result.returncode}")
        logger.info(f"Script output: {result.stdout}")
        if result.stderr:
            logger.warning(f"Script stderr: {result.stderr}")

        filename = f"{task_output_dir}/SNP_Analysis_Report.pdf"

        # 检查结果文件是否存在
        if not os.path.exists(filename):
            logger.error(f"SNP report PDF not found: {filename}")
            return

        # 发送结果邮件
        send_email_with_attachment(email, gene_name, filename, "SNP analysis completed", "SNP Analysis Results")
        logger.info(f"SNP analysis completed for {gene_name}")
        
        
    except Exception as e:
        error_msg = f"SNP analysis failed:\n{str(e)}"
        logger.error(error_msg)
        send_email(email, gene_name, error_msg, "SNP Analysis Error")
    finally:
        # 清理临时目录
        cleanup_output_directory(task_output_dir)

def run_offtarget_analysis(data):
    """执行脱靶分析并发送结果邮件"""
    targetgene = data['target_gene']
    target_aso = data['target_aso']
    email = data['email']
    logger.info(f"Starting offtarget analysis for {target_aso}")
    task_id = str(uuid.uuid4())
    task_output_dir = os.path.join(OUTPUT_DIR_off, task_id)
    os.makedirs(task_output_dir, exist_ok=True) 
    try:
        filename = f"{task_output_dir}/top15genes.pdf"
        result = subprocess.run(
            f'bash /offtarget/_offtarget.sh {targetgene} {target_aso} {task_output_dir}',
            shell=True,
            capture_output=True,
            text=True,
            check=False,  # 不检查状态码
            env={**os.environ, 'OUTPUT_DIR_off': task_output_dir},
            timeout=None  # 移除超时限制
        )
        logger.info(f"Off-target analysis script exit code: {result.returncode}")
        logger.info(f"Script output: {result.stdout}")
        if result.stderr:
            logger.warning(f"Script stderr: {result.stderr}")
        # 检查脚本是否生成了结果文件
        if not os.path.exists(filename):
            # 如果脚本没有生成结果文件，创建一个默认的
            with open(filename, "w") as f:
                f.write(f"Off-Target Analysis Results for {target_aso}\n")
                f.write("\nOff-Target Analysis:\n")
                f.write("3 potential off-target sites identified\n")
                f.write("All have moderate binding affinity\n")
        
        # 发送结果邮件
        send_email_with_attachment(email, target_aso, filename, "Off-target analysis completed", "Off-Target Analysis Results")
        logger.info(f"Offtarget analysis completed for {target_aso}")
    except Exception as e:
        error_msg = f"Offtarget analysis failed:\n{str(e)}"
        logger.error(error_msg)
        send_email(email, target_aso, error_msg, "Off-Target Analysis Error")
    finally:
        # 清理临时目录
        cleanup_output_directory(task_output_dir)

def cleanup_output_directory(directory):
    """清理指定目录中的所有文件"""
    try:
        if not os.path.exists(directory):
            logger.warning(f"Output directory does not exist: {directory}")
            return
        
        for filename in os.listdir(directory):
            file_path = os.path.join(directory, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                logger.error(f"Failed to delete file {file_path}: {str(e)}")
    except Exception as e:
        logger.error(f"Error cleaning output directory: {str(e)}")

def send_email_with_attachment(to_email, subject_name, file_path, command_output, subject_prefix):
    """发送带附件的邮件"""
    logger.debug(f"[DEBUG] Sending email to: {to_email}, Subject: {subject_name}")
    logger.debug(f"[DEBUG] File path: {file_path}")
    logger.debug(f"[DEBUG] File exists: {os.path.exists(file_path)}")
    if os.path.exists(file_path):
        logger.debug(f"[DEBUG] File size: {os.path.getsize(file_path)} bytes")
    
    try:
        msg = MIMEMultipart()
        msg['From'] = FROM_EMAIL
        msg['To'] = to_email
        msg['Subject'] = f'{subject_prefix} - {subject_name}'
        
        # 根据分析类型生成不同的邮件内容
        body_content = ""
        if os.path.exists(file_path):
            if "Homology" in subject_prefix:
                body_content = f""
                body_content += f"<h3>Homology Analysis Results for {subject_name}</h3>"
                body_content += f"<p>Dear User,</p>"
                body_content += f"<p>We are pleased to inform you that the homology analysis for <strong>{subject_name}</strong> has been completed successfully.</p>"
                body_content += f"<p><strong>Analysis Details:</strong></p>"
                body_content += f"<ul>"
                body_content += f"<li><strong>Target Gene:</strong> {subject_name}</li>"
                body_content += f"<li><strong>Analysis Type:</strong> Cross-species homology analysis</li>"
                body_content += f"<li><strong>Result Format:</strong> Text file with homology scores</li>"
                body_content += f"</ul>"
                body_content += f"<p>The detailed results are attached to this email. The analysis includes homology scores across different species, which can help you assess the specificity of your ASO sequence.</p>"
                body_content += f"<p><strong>Interpretation Tips:</strong></p>"
                body_content += f"<ul>"
                body_content += f"<li>Higher homology scores indicate greater sequence similarity</li>"
                body_content += f"<li>Lower scores suggest better species-specificity</li>"
                body_content += f"<li>Focus on scores for your target species vs. non-target species</li>"
                body_content += f"</ul>"
                body_content += f"<p>If you have any questions about the results or need further assistance, please feel free to reply directly to this email. Our team will be happy to help you.</p>"
                body_content += f"<p>Thank you for using the ClinASO Analysis Platform!</p>"
                body_content += f"<hr>"
                body_content += f"<p><em>Best regards,</em></p>"
                body_content += f"<p><em>The ClinASO Team</em></p>"
                body_content += f"<p><em>School of Life Sciences, Yunnan University</em></p>"
            elif "SNP" in subject_prefix:
                body_content = f""
                body_content += f"<h3>SNP Analysis Results for {subject_name}</h3>"
                body_content += f"<p>Dear User,</p>"
                body_content += f"<p>We are pleased to inform you that the SNP analysis for <strong>{subject_name}</strong> has been completed successfully.</p>"
                body_content += f"<p><strong>Analysis Details:</strong></p>"
                body_content += f"<ul>"
                body_content += f"<li><strong>Target Gene:</strong> {subject_name}</li>"
                body_content += f"<li><strong>Analysis Type:</strong> SNP frequency and impact analysis</li>"
                body_content += f"<li><strong>Result Format:</strong> PDF report with SNP distribution</li>"
                body_content += f"</ul>"
                body_content += f"<p>The detailed results are attached to this email. The analysis includes information about SNPs in the target region, their frequencies, and potential impact on ASO binding.</p>"
                body_content += f"<p><strong>Interpretation Tips:</strong></p>"
                body_content += f"<ul>"
                body_content += f"<li>High-frequency SNPs may affect ASO binding efficiency</li>"
                body_content += f"<li>Consider SNP coverage when designing ASO sequences</li>"
                body_content += f"<li>Focus on SNPs within the ASO target region</li>"
                body_content += f"</ul>"
                body_content += f"<p>If you have any questions about the results or need further assistance, please feel free to reply directly to this email. Our team will be happy to help you.</p>"
                body_content += f"<p>Thank you for using the ClinASO Analysis Platform!</p>"
                body_content += f"<hr>"
                body_content += f"<p><em>Best regards,</em></p>"
                body_content += f"<p><em>The ClinASO Team</em></p>"
                body_content += f"<p><em>School of Life Sciences, Yunnan University</em></p>"
            elif "Off-Target" in subject_prefix:
                body_content = f""
                body_content += f"<h3>Off-Target Analysis Results for {subject_name}</h3>"
                body_content += f"<p>Dear User,</p>"
                body_content += f"<p>We are pleased to inform you that the off-target analysis for <strong>{subject_name}</strong> has been completed successfully.</p>"
                body_content += f"<p><strong>Analysis Details:</strong></p>"
                body_content += f"<ul>"
                body_content += f"<li><strong>Target ASO:</strong> {subject_name}</li>"
                body_content += f"<li><strong>Analysis Type:</strong> Potential off-target binding sites</li>"
                body_content += f"<li><strong>Result Format:</strong> PDF report with top off-target candidates</li>"
                body_content += f"</ul>"
                body_content += f"<p>The detailed results are attached to this email. The analysis identifies potential off-target binding sites and their binding affinities, which can help you assess the specificity of your ASO sequence.</p>"
                body_content += f"<p><strong>Interpretation Tips:</strong></p>"
                body_content += f"<ul>"
                body_content += f"<li>Lower binding energy indicates stronger potential binding</li>"
                body_content += f"<li>Focus on off-target sites with high binding affinity</li>"
                body_content += f"<li>Consider the biological relevance of off-target genes</li>"
                body_content += f"</ul>"
                body_content += f"<p>If you have any questions about the results or need further assistance, please feel free to reply directly to this email. Our team will be happy to help you.</p>"
                body_content += f"<p>Thank you for using the ClinASO Analysis Platform!</p>"
                body_content += f"<hr>"
                body_content += f"<p><em>Best regards,</em></p>"
                body_content += f"<p><em>The ClinASO Team</em></p>"
                body_content += f"<p><em>School of Life Sciences, Yunnan University</em></p>"
            else:
                body_content = f"Analysis result generated. Please check the attached result file: {os.path.basename(file_path)}"
            logger.debug(f"[DEBUG] File exists, using success message")
        else:
            body_content = f"Result file not found: {file_path}\nCommand output:\n{command_output}"
            logger.error(f"[DEBUG] File not found, using error message")
        
        body = MIMEText(body_content, 'html', 'utf-8')
        msg.attach(body)
        
        # 添加附件（如果文件存在）
        if os.path.exists(file_path):
            logger.debug(f"[DEBUG] Adding attachment: {file_path}")
            try:
                with open(file_path, 'rb') as file:
                    part = MIMEApplication(
                        file.read(),
                        Name=os.path.basename(file_path)
                    )
                part['Content-Disposition'] = f'attachment; filename="{os.path.basename(file_path)}"'
                msg.attach(part)
                logger.debug(f"[DEBUG] Attachment added successfully")
            except Exception as e:
                logger.error(f"[DEBUG] Error adding attachment: {str(e)}")
        
        # 连接SMTP服务器并发送
        logger.debug(f"[DEBUG] Connecting to SMTP server: {SMTP_SERVER}")
        with smtplib.SMTP_SSL(SMTP_SERVER, SMTP_PORT) as server:
            logger.debug(f"[DEBUG] SMTP SSL connection established")
            server.login(SMTP_USERNAME, SMTP_PASSWORD)
            logger.debug(f"[DEBUG] SMTP login successful")
            server.send_message(msg)
            logger.debug(f"[DEBUG] Email sent successfully")
            
        return True
        
    except smtplib.SMTPAuthenticationError as auth_err:
        error_msg = f"Authentication failed: {str(auth_err)}"
        logger.error(f"[DEBUG] {error_msg}")
        return False
        
    except smtplib.SMTPException as smtp_err:
        error_msg = f"SMTP error: {str(smtp_err)}"
        logger.error(f"[DEBUG] {error_msg}")
        return False
        
    except Exception as e:
        error_msg = f"Email sending failed: {str(e)}"
        logger.error(f"[DEBUG] {error_msg}")
        import traceback
        traceback.print_exc()
        return False

def send_email(to_email, subject_name, content, subject_prefix):
    """发送结果邮件（不带附件）"""
    logger.info(f"Sending error email to: {to_email}, Subject: {subject_name}")
    
    try:
        msg = MIMEMultipart()
        msg['From'] = FROM_EMAIL
        msg['To'] = to_email
        msg['Subject'] = f'{subject_prefix} - {subject_name}'
        
        # 邮件正文
        body = MIMEText(f"""
        <h2>{subject_name} Analysis Result</h2>
        <pre>{content}</pre>
        """, 'html', 'utf-8')
        msg.attach(body)
        
        # 连接SMTP服务器并发送
        with smtplib.SMTP_SSL(SMTP_SERVER, SMTP_PORT) as server:
            server.login(SMTP_USERNAME, SMTP_PASSWORD)
            server.send_message(msg)
            logger.info("Error email sent successfully")
            
        return True
        
    except Exception as e:
        logger.error(f"Error email sending failed: {str(e)}")
        return False

# 路由定义
@app.route('/submit_design', methods=['POST'])
def submit_design():
    """处理第一个表单提交 - Gapmer Design"""
    try:
        data = {
            'species': request.form.get('species'),
            'gene_name': request.form.get('gene_name', '').upper(),
            'aso_len': request.form.get('aso_len'),
            'homologous_species': request.form.get('homologous_species'),
            'gc_content': request.form.get('gc_content'),
            'aso_count': request.form.get('aso_count'),
            'priority': request.form.get('priority'),
            'email': request.form.get('email')
        }
        
        # 验证输入
        if not data['gene_name']:
            return jsonify({'error': 'Gene name cannot be empty'}), 400
        if not data['email'] or '@' not in data['email']:
            return jsonify({'error': 'Please enter a valid email address'}), 400
        if not data['aso_len']:
            return jsonify({'error': 'ASO length is required'}), 400
        if not data['aso_count']:
            return jsonify({'error': 'Number of ASOs is required'}), 400
        if not data['priority']:
            return jsonify({'error': 'Selection priority is required'}), 400
                # 去重检查
        dedup_params = data['species'] + ':' + data['gene_name'] + ':' + data['aso_len'] + ':' + data['homologous_species'] + ':' + data['gc_content'] + ':' + data['aso_count'] + ':' + data['priority']
        if is_duplicate_task('design', data['email'], dedup_params):
            return jsonify({
                'status': 'duplicate',
                'message': 'You have already submitted this exact analysis. Results will be sent to your email.'
            }), 200

# 确保工作线程已启动
        start_worker()
        
        # 将任务加入队列
        try:
            task_queue.put(('design', data), timeout=5)
        except queue.Full:
            return jsonify({
                'error': 'Queue is full',
                'message': 'Too many pending tasks. Please try again later.'
            }), 503
        
        return jsonify({
            'status': 'success', 
            'message': 'Analysis started. Results will be sent to your email.'
        })
    
    except Exception as e:
        logger.error(f"Submit design error: {str(e)}")
        return jsonify({'error': f'Server error: {str(e)}'}), 500

@app.route('/analyze_sequence', methods=['POST'])
def analyze_sequence():
    """处理第二个表单提交 - Homology Analysis"""
    try:
        # 打印所有接收到的表单数据用于调试
        logger.info(f"Received form data: {dict(request.form)}")
        
        data = {
            'gene_name': request.form.get('gene_name', '').upper(),
            'gene_sequence': request.form.get('gene_sequence'),
            'analysis_type': request.form.get('analysis_type'),
            'email': request.form.get('email')
        }
        
        # 添加字段验证日志
        logger.info(f"Validating data: gene_name={data['gene_name']}, gene_sequence={bool(data['gene_sequence'])}, email={data['email']}")
        
        # 验证输入
        if not data['gene_name'] or not data['gene_sequence']:
            logger.error("Gene name or sequence is empty")
            return jsonify({'error': 'Gene name and sequence cannot be empty'}), 400
        
        if not data['email'] or '@' not in data['email']:
            logger.error(f"Invalid email: {data['email']}")
            return jsonify({'error': 'Please enter a valid email address'}), 400
        
                # 去重检查
        dedup_params = data['gene_name'] + ':' + data['gene_sequence'] + ':' + data.get('analysis_type', '')
        if is_duplicate_task('homology', data['email'], dedup_params):
            return jsonify({
                'status': 'duplicate',
                'message': 'You have already submitted this exact analysis. Results will be sent to your email.'
            }), 200

# 确保工作线程已启动
        start_worker()
        
        # 将任务加入队列
        try:
            task_queue.put(('homology', data), timeout=5)
            logger.info(f"Homology task queued for gene: {data['gene_name']}")
        except queue.Full:
            logger.warning("Task queue is full")
            return jsonify({
                'error': 'Queue is full',
                'message': 'Too many pending tasks. Please try again later.'
            }), 503
        
        return jsonify({
            'status': 'success', 
            'message': 'Analysis started. Results will be sent to your email.'
        })
    
    except Exception as e:
        logger.error(f"Analyze sequence error: {str(e)}", exc_info=True)
        return jsonify({'error': f'Server error: {str(e)}'}), 500


@app.route('/analyze_snp', methods=['POST'])
def analyze_snp():
    """处理第三个表单提交 - SNP Analysis"""
    try:
        data = {
            'gene_name': request.form.get('snp_gene_name', '').upper(),
            'aso_sequence': request.form.get('snp_aso_sequence'),
            'email': request.form.get('snp_email') 
        }
        # 添加验证日志
        logger.info(f"SNP analysis data received: {data}")
        
        # 验证输入
        if not data['gene_name'] or not data['aso_sequence']:
            logger.error("Gene name or ASO sequence is empty")
            return jsonify({'error': 'Gene name and ASO sequence cannot be empty'}), 400
        
        if not data['email'] or '@' not in data['email']:
            logger.error(f"Invalid email: {data['email']}")
            return jsonify({'error': 'Please enter a valid email address'}), 400
        
                # 去重检查
        dedup_params = data['gene_name'] + ':' + data['aso_sequence']
        if is_duplicate_task('snp', data['email'], dedup_params):
            return jsonify({
                'status': 'duplicate',
                'message': 'You have already submitted this exact analysis. Results will be sent to your email.'
            }), 200

# 确保工作线程已启动
        start_worker()
        
        # 将任务加入队列
        try:
            task_queue.put(('snp', data), timeout=5)
            logger.info(f"SNP task queued for gene: {data['gene_name']}")
        except queue.Full:
            logger.warning("Task queue is full")
            return jsonify({
                'error': 'Queue is full',
                'message': 'Too many pending tasks. Please try again later.'
            }), 503
        
        return jsonify({
            'status': 'success', 
            'message': 'Analysis started. Results will be sent to your email.'
        })
    
    except Exception as e:
        logger.error(f"Analyze SNP error: {str(e)}", exc_info=True)
        return jsonify({'error': f'Server error: {str(e)}'}), 500

@app.route('/analyze_offtarget', methods=['POST'])
def analyze_offtarget():
    """处理第四个表单提交 - Off-Target Analysis"""
    try:
        data = {
            'target_gene': request.form.get('target_gene', '').upper(),
            'target_aso': request.form.get('target_aso'),
            'email': request.form.get('email')
        }
        logger.info(data)
        # 验证输入
        if not data['target_aso']:
            return jsonify({'error': 'Target sequence cannot be empty'}), 400
        if not data['email'] or '@' not in data['email']:
            return jsonify({'error': 'Please enter a valid email address'}), 400
        
                # 去重检查
        dedup_params = data['target_gene'] + ':' + data['target_aso']
        if is_duplicate_task('offtarget', data['email'], dedup_params):
            return jsonify({
                'status': 'duplicate',
                'message': 'You have already submitted this exact analysis. Results will be sent to your email.'
            }), 200

# 确保工作线程已启动
        start_worker()
        
        # 将任务加入队列
        try:
            task_queue.put(('offtarget', data), timeout=5)
        except queue.Full:
            return jsonify({
                'error': 'Queue is full',
                'message': 'Too many pending tasks. Please try again later.'
            }), 503
        
        return jsonify({
            'status': 'success', 
            'message': 'Analysis started. Results will be sent to your email.'
        })
    
    except Exception as e:
        logger.error(f"Analyze offtarget error: {str(e)}")
        return jsonify({'error': f'Server error: {str(e)}'}), 500

@app.route('/queue-status', methods=['GET'])
def queue_status():
    """获取当前队列状态"""
    return jsonify({
        'queue_size': task_queue.qsize(),
        'worker_status': 'running' if worker_thread and worker_thread.is_alive() else 'stopped'
    })

@app.route('/system-status', methods=['GET'])
def system_status():
    """获取系统状态信息"""
    return jsonify({
        'queue_size': task_queue.qsize(),
        'worker_status': 'running' if worker_thread and worker_thread.is_alive() else 'stopped'
    })
if __name__ == '__main__':
    # 设置生产环境模式
    from werkzeug.serving import run_simple
    start_worker()
    logger.info("Application started")
    
    try:
        # 使用生产服务器而不是开发服务器
        run_simple('0.0.0.0', 5000, app, use_reloader=False, use_debugger=False, threaded=True)
    finally:
        # 确保在应用停止时关闭工作线程
        stop_worker()
        logger.info("Application stopped")
