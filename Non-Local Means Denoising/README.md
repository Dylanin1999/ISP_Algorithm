1、Denoise.cpp: 对Non-Local Means的实现

​	 ①void NonLocalMeans():  没有进行优化的NLM算法，时间复杂度高，速度慢。对Non-local Means算法的原始想法进行了最基础还原，想要学习和了解Non-local Means算法的,可以直接查看该文件。

​	②void NLMWithIntegralImg():　使用积分图对计算步骤进行加速，速度极大地提升。

2、IntegralImg.h：对积分图的计算

​	利用了特殊方法对积分图进行计算，计算过程更为简洁明了。

