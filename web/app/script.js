// 选项卡切换功能
document.querySelectorAll('.tab-btn').forEach(button => {
    button.addEventListener('click', () => {
        // 移除所有活动状态
        document.querySelectorAll('.tab-btn').forEach(btn => {
            btn.classList.remove('active');
        });
        document.querySelectorAll('.tab-content').forEach(content => {
            content.classList.remove('active');
        });
        
        // 添加当前活动状态
        button.classList.add('active');
        const tabId = button.getAttribute('data-tab');
        document.getElementById(tabId).classList.add('active');
    });
});

// 表单提交处理
document.addEventListener('DOMContentLoaded', function() {
    // 表单1: Gapmer Design (邮件发送结果)
    const designForm = document.getElementById('analysis-form');
    if (designForm) {
        designForm.addEventListener('submit', function(e) {
            e.preventDefault();
            handleFormSubmit(this, '#response-message', '/submit_design');
        });
    }
    
    // 表单2: Homology Analysis (邮件发送结果) 
    const sequenceForm = document.getElementById('sequence-form');
    if (sequenceForm) {
        sequenceForm.addEventListener('submit', function(e) {
            e.preventDefault();
            handleFormSubmit(this, '#sequence-response', '/analyze_sequence');
        });
    }

    // 表单3: SNP Analysis (邮件发送结果)
    const snpForm = document.getElementById('snp-form');
    if (snpForm) {
        snpForm.addEventListener('submit', function(e) {
            e.preventDefault();
            handleFormSubmit(this, '#snp-response', '/analyze_snp');
        });
    }
    
    // 表单4: Off-Target Analysis (邮件发送结果)
    const offtargetForm = document.getElementById('offtarget-form');
    if (offtargetForm) {
        offtargetForm.addEventListener('submit', function(e) {
            e.preventDefault();
            handleFormSubmit(this, '#offtarget-response', '/analyze_offtarget');
        });
    }
    
    // DNA动画效果
    const dnaElement = document.querySelector('.dna-animation');
    if (dnaElement) {
        const symbols = ['🧬', '🔬', '🧪', '🧫'];
        let currentIndex = 0;
        setInterval(() => {
            dnaElement.textContent = symbols[currentIndex];
            currentIndex = (currentIndex + 1) % symbols.length;
        }, 2500);
    }
});

// GC内容验证函数
function validateGCContent(input) {
    const value = parseFloat(input.value);
    if (value < 0.35 || value > 0.7) {
        input.setCustomValidity('GC content must be between 0.35 and 0.7');
    } else {
        input.setCustomValidity('');
    }
}

// 通用表单提交处理函数（所有表单都通过邮件发送结果）
async function handleFormSubmit(form, responseSelector, endpoint) {
    const formData = new FormData(form);
    const submitBtn = form.querySelector('.btn-submit');
    const responseMessage = document.querySelector(responseSelector);
    
    // 显示加载状态
    const originalText = submitBtn.textContent;
    submitBtn.textContent = 'Processing...';
    submitBtn.disabled = true;
    responseMessage.style.display = 'none';
    
    try {
        console.log('Submitting form to:', endpoint);
        console.log('Form data:', Object.fromEntries(formData));
        
        const response = await fetch(endpoint, {
            method: 'POST',
            body: formData
        });
        
        console.log('Response status:', response.status);
        console.log('Response status text:', response.statusText);
        console.log('Response headers:', response.headers);
        
        // 检查内容类型
        const contentType = response.headers.get('content-type');
        console.log('Content-Type:', contentType);
        
        // 只读取一次响应体
        const responseText = await response.text();
        console.log('Response text:', responseText);
        
        // 尝试解析JSON
        let responseData;
        if (contentType && contentType.includes('application/json')) {
            try {
                responseData = JSON.parse(responseText);
            } catch (e) {
                console.error('Failed to parse JSON response:', e);
                throw new Error('Invalid JSON format in server response');
            }
        } else {
            // 如果不是JSON，创建一个简单的响应对象
            responseData = {
                message: responseText,
                error: response.ok ? null : 'Server returned non-JSON response'
            };
            
            // 检查是否是HTML错误页面
            if (responseText.trim().startsWith('<!DOCTYPE') || responseText.trim().startsWith('<html')) {
                console.warn('Server returned HTML response instead of JSON');
                // 尝试提取错误信息
                const titleMatch = responseText.match(/<title>(.*?)<\/title>/i);
                const title = titleMatch ? titleMatch[1] : 'Unknown Error';
                responseData.error = `Server error: ${title}`;
                responseData.message = 'The server encountered an error. Please try again later.';
            }
        }
        
        if (!response.ok) {
            throw new Error(responseData.error || `HTTP error! status: ${response.status}`);
        }
        
        console.log('Response data:', responseData);
        
        if (responseData.status === 'duplicate') {
            responseMessage.textContent = '⚠️ ' + (responseData.message || 'You have already submitted this analysis. Results will be sent to your email.');
            responseMessage.className = 'response-message warning';
        } else if (responseData.error) {
            responseMessage.textContent = responseData.error + (responseData.message ? ': ' + responseData.message : '');
            responseMessage.className = 'response-message error';
        } else {
            responseMessage.textContent = responseData.message;
            responseMessage.className = 'response-message success';
        }
        
        responseMessage.style.display = 'block';
        
        // 重置表单（除了第一个表单，因为它包含email字段）
        if (form.id !== 'analysis-form') {
            form.reset();
        }
    } catch (error) {
        console.error('Form submission error:', error);
        responseMessage.textContent = 'Error: ' + error.message;
        responseMessage.className = 'response-message error';
        responseMessage.style.display = 'block';
    } finally {
        submitBtn.textContent = originalText;
        submitBtn.disabled = false;
    }
}

// ASO count validation
function validateAsoCount(input) {
    const value = parseInt(input.value);
    if (value > 1500) {
        input.setCustomValidity('Number of ASOs cannot exceed 1500');
    } else if (value < 1) {
        input.setCustomValidity('Number of ASOs must be at least 1');
    } else {
        input.setCustomValidity('');
    }
}

