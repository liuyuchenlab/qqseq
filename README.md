# qqseq！ 
###### 在学习基因敲除实验时写的，可以一键获得转录本的外显子和内含子序列


library(devtools)  


###### 如果连接失败：  
###### 1.尝试修复Hosts配置  
###### 2.尝试Win+R，输入inetcpl.cpl 直接打开Internet选项。打开后，在高级中勾选使用TLS 1.0、使用TLS 1.1、使用TLS 1.2、使用TLS 1.3。


devtools::install_github('liuyuchenlab/qqseq')  
library(qqseq) 


# 示例使用
###### 使用esembl的转录本id来进行索引，目前只能用mm39的小鼠和hg38的人类，需要其他物种或者参考基因组的可以自己改一下代码，十分方便
###### 小鼠mouse，人类human


result <- qqseq("ENSMUST00000127786", "mouse")


###### 这里会在当前文件夹保存一个transcrpt_id.csv的文件
#### 结果与ESEMBL数据库一致


