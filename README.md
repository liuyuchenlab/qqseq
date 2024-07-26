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
#### start和end与ESEMBL是反的
#### 保存的序列与ESEMBL一致


# 负义链
##### ENSMUST00000153883
###### ESEMBL
![image](https://github.com/user-attachments/assets/c76753ae-77bd-4c38-9ce4-8357c6bdc872)


###### qqseq位置
![image](https://github.com/user-attachments/assets/f88d7fe3-5ff9-4e93-88bb-c69ddd40d819)


###### qqseq序列
![image](https://github.com/user-attachments/assets/5d917bce-7129-4b3f-aea4-553627fa1e4e)



# 正义链
##### ENSMUST00000152916
###### ESEMBL
![image](https://github.com/user-attachments/assets/5533702a-b1b0-485c-a4a2-24373c6a64e3)


###### qqseq位置
![image](https://github.com/user-attachments/assets/602a0a03-f713-4ef9-84ed-86fe07f4f21a)


###### qqseq序列
![image](https://github.com/user-attachments/assets/52d9af0b-53b8-4cab-9a1d-6386be75816b)



###### 有时会出现esembl链接超时的问题，不过问题不大，再试一次就好了



