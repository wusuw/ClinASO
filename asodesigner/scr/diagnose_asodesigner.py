#!/usr/bin/env python3
"""
ASO Designer 后端诊断工具
确保在拥有sudo权限的用户下运行
"""

import os
import subprocess
import logging
from datetime import datetime

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    filename=f"asodesigner_diagnose_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
)
logger = logging.getLogger("ASO_Designer_Diagnose")

def run_command(cmd, description):
    """执行命令并记录结果"""
    logger.info(f"执行: {cmd}")
    try:
        result = subprocess.run(
            cmd, 
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        logger.info(f"✅ {description}成功")
        logger.debug(f"输出:\n{result.stdout}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"❌ {description}失败")
        logger.error(f"错误: {e.stderr}")
        return False

def main():
    # 1. 检查关键目录权限
    required_dirs = [
        "/asodesigner/scr",
        "/asodesigner/scr/results",
        "/asodesigner/scr/logs"
    ]
    
    for directory in required_dirs:
        if not os.path.exists(directory):
            run_command(f"sudo mkdir -p {directory}", f"创建目录 {directory}")
        
        run_command(f"sudo ls -ld {directory}", f"检查目录权限 {directory}")
        run_command(f"sudo chmod 775 {directory}", f"设置目录权限 {directory}")
    
    # 2. 检查脚本位置与权限
    script_path = "/asodesigner/scr/aso_design.sh"
    
    # 模拟创建脚本（如果不存在）
    if not os.path.exists(script_path):
        logger.warning("aso_design.sh不存在，创建模拟脚本")
        with open(script_path, "w") as f:
            f.write("#!/bin/bash\necho '模拟ASO设计脚本运行中...'\nsleep 3\necho '任务完成!' > results.txt")
        
        run_command(f"sudo chmod +x {script_path}", "设置脚本权限")
    
    run_command(f"ls -l {script_path}", "检查脚本权限")
    
    # 3. 测试环境依赖
    run_command("python3 --version", "检查Python版本")
    run_command("pip3 list | grep Flask", "检查Flask安装")
    run_command("which gunicorn || echo 'Gunicorn未安装'", "检查Gunicorn")
    
    # 4. 测试服务启动
    print("\n尝试启动后端服务...")
    os.chdir("/asodesigner/scr")
    run_command("sudo python3 app.py", "启动后端服务")

if __name__ == "__main__":
    main()
    print("\n诊断完成！详细日志已保存到当前目录")
