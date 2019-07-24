# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 08:03:54 2019

@author: LJ
"""

#! /bin/bash
usage() { echo "Usage: $0 [-r <电话号码>] [-c <短信内容>] [-b <营销短信请设置-b参数>]" 1>&2; exit 1; }
APP_ID="free_trial"
TOKEN="MjYxOTNjMTkyZTJjZjgzODA5OGVkNjQyYzgzOGUwMjY="
USER_ID="0"
RECEIVER=""
DESCRIPTION=""
BUSINESS=false
while getopts 'r:c:b' OPT;
do
   case $OPT in
     r)
        RECEIVER="$OPTARG";;
     c)
        DESCRIPTION="$OPTARG";;
     b)
        BUSINESS=true;;
     *)
        usage
   esac
done
shift $(($OPTIND - 1))
if [ -z "$RECEIVER" ] || [ -z "$DESCRIPTION" ] || [ "-b" == "$DESCRIPTION" ]; then
    usage
fi
DATA="[{channel : \"sms\", description : \"$DESCRIPTION\",receiver : \"$RECEIVER\",business : \"$BUSINESS\"}]"
INPUT=${TOKEN}${APP_ID}${DATA}
SIGNATURE=`echo -n $INPUT|md5sum|awk '{print $1}'`
curl -H "appid:$APP_ID" -H "token:$TOKEN" -H "userid:$USER_ID" -H "signature:$SIGNATURE" -X POST -d "$DATA" http://gaojing.baidu.com/AlertList/push
#--------------------- 
#作者：百度告警 
#来源：CSDN 
#原文：https://blog.csdn.net/baidu_tonggao/article/details/50429201 
#版权声明：本文为博主原创文章，转载请附上博文链接！