# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:56:15 2019

@author: 我永远爱八重樱
"""

import smtplib
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart

def send_mail(ipaddr,news=''):
    msg_from='635080302@qq.com'                                 #发送方邮箱
    passwd='kpkeottafbirbead'                                   #填入发送方邮箱的授权码
    msg_to='635080302@qq.com'                                  #收件人邮箱
                                
    subject=ipaddr
    content=news # 正文
    msg = MIMEMultipart('related')
    m_content = MIMEText(content,'html','utf-8')
    msg.attach(m_content)
    msg['Subject'] = subject
    msg['From'] = msg_from
    msg['To'] = msg_to
    try:
        s = smtplib.SMTP_SSL("smtp.qq.com",465) #邮件服务器及端口号
        s.login(msg_from, passwd)
        s.sendmail(msg_from, msg_to, msg.as_string())
        print("发送成功")
    except:
        print("发送失败")
    finally:
        s.quit()

send_mail('test')