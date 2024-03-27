import os
import cv2   
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import itertools
# with open("x09.raw","rb") as f:
#     data= f.read()  
#     print(len(data))
#     data = list(bytes(data))
#     data = np.array(data)
#     data = data.astype(np.float16)
    
#     # data = np.reshape(3848,2168)
#     data = data.reshape(3848*2176,2)
   
#     img = (data[:,1])*256+data[:,0]

#     img = np.array(img)
#     img = img.reshape(2176,3848)
    

#     # im = Image.fromarray(img)

#     # print(img)
#     # print("=======================")

#     # #转成低位后保存
#     # img = img/256
#     # cv2.imwrite('img.jpeg', img)


#     #奇数行
#     odd_rows = np.array([row for index, row in enumerate(img) if index % 2 == 0])
#     #偶数行
#     even_rows = np.array([row for index, row in enumerate(img) if index % 2 != 0])

#     #提取R分量 GR分量
#     transpose_rows = odd_rows

#     transpose_rows = np.array(np.transpose(odd_rows))
#     print(transpose_rows.shape)
#     print(odd_rows.shape)
#     R_Vec = np.array([row for index, row in enumerate(transpose_rows) if index % 2 == 0])
#     GR_Vec = np.array([row for index, row in enumerate(transpose_rows) if index % 2 != 0])
#     R_Vec = np.array(np.transpose(R_Vec))
#     GR_Vec = np.array(np.transpose(GR_Vec))
#     print(R_Vec.shape)





#     # odd_rows_img = odd_rows/256
#     # odd_rows_img=odd_rows_img.astype(np.uint8)   #进行类型转换

#     # Rvec_rows_img = R_Vec/256
#     # Rvec_rows_img=Rvec_rows_img.astype(np.uint8)   #进行类型转换

#     # cv2.imshow("odd",Rvec_rows_img)
#     # cv2.imwrite('odd.jpeg', odd_rows/256)


#     #提取B分量 GB分量
#     transpose_rows = even_rows

#     transpose_rows = np.array(np.transpose(transpose_rows))
#     print(transpose_rows.shape)
#     print(odd_rows.shape)
#     GB_Vec = np.array([row for index, row in enumerate(transpose_rows) if index % 2 == 0])
#     B_Vec = np.array([row for index, row in enumerate(transpose_rows) if index % 2 != 0])
#     GB_Vec = np.array(np.transpose(GB_Vec))
#     B_Vec = np.array(np.transpose(B_Vec))

#     Rvec_rows_img = R_Vec/256
#     Bvec_rows_img = B_Vec/256
#     Gvec_rows_img = (GB_Vec+GR_Vec)/2/256 

#     Rvec_rows_img=Rvec_rows_img.astype(np.uint8)   #进行类型转换
#     Gvec_rows_img=Gvec_rows_img.astype(np.uint8)   #进行类型转换
#     Bvec_rows_img=Bvec_rows_img.astype(np.uint8)   #进行类型转换

#     merged = cv2.merge([Bvec_rows_img,Gvec_rows_img,Rvec_rows_img])#合并R、G、B分量
    
#     # img = merged/256

#     print("imgsize is",img.shape)

#     # cv2.imshow("Merged",img)
#     # cv2.imshow("odd",R_Vec/256)
#     cv2.imshow('Rvec_rows_img.jpeg', Rvec_rows_img)
#     cv2.imshow('Gvec_rows_img.jpeg', Gvec_rows_img)
#     cv2.imshow('Bvec_rows_img.jpeg', Bvec_rows_img)
#     cv2.imshow('odd.jpeg', merged)

#     # print(odd_rows)
#     # cv2.imwrite('odd.jpeg', merged)
#     cv2.waitKey(0)
#     # print("=======================")

#     # print(even_rows)

#     # data = data.reshape(3840,2160)
#     # print(data)


img = [
        [1,2,1,2,1,2,1,2],
        [3,4,3,4,3,4,3,4],
        [1,2,1,2,1,2,1,2],
        [3,4,3,4,3,4,3,4]
        ]

img_1 = np.zeros(4*4).reshape(4,4)
padding_1 = np.zeros(4*4).reshape(4,4)

padding_2 = np.zeros(4*8).reshape(4,8)

img_1 = img_1+1
print(img_1)
# img_1 = itertools.chain.from_iterable(zip(padding_1,img_1))
img_1 = np.array(list(itertools.chain.from_iterable(zip(img_1,padding_1))))
print(np.array(list(img_1)))
print(img_1.shape)
# img_1 = itertools.chain.from_iterable(zip(img_1,padding))
#完成重新排布
img_1 = np.array(list(itertools.chain.from_iterable(zip(np.array(np.transpose(img_1)),padding_2))))
print(np.array(list(img_1)))

img = np.array(img)
def extract_rgb(img,cfa_pattern):

    size_height = int(img.shape[0])
    size_width = int(img.shape[1])
    print(size_height,size_width)
    CFA_0_0 = np.zeros(size_height*size_width).reshape(size_height,size_width)
    CFA_0_1 = np.zeros(size_height*size_width).reshape(size_height,size_width)
    CFA_1_0 = np.zeros(size_height*size_width).reshape(size_height,size_width)
    CFA_1_1 = np.zeros(size_height*size_width).reshape(size_height,size_width)

    #split rggb
    for vertical_idx in range(0,size_height//2):
        for horizon_idx in range(0,size_width//2):
            CFA_0_0[0+vertical_idx*2][horizon_idx*2] = img[0+vertical_idx*2][horizon_idx*2]
            CFA_0_1[0+vertical_idx*2][horizon_idx*2+1] = img[0+vertical_idx*2][horizon_idx*2+1]
            CFA_1_0[0+vertical_idx*2+1][horizon_idx*2] = img[0+vertical_idx*2+1][horizon_idx*2]
            CFA_1_1[0+vertical_idx*2+1][horizon_idx*2+1] = img[0+vertical_idx*2+1][horizon_idx*2+1]





#     # if cfa_pattern=="rggb":


extract_rgb(img,"rggb")