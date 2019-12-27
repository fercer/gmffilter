import cv2
import numpy as np
import matplotlib.pyplot as plt
import gmf

import os
import sys
import time

print('Testing GMF filter')

img_cat = np.zeros([0, 300, 300])
img_gt_cat = np.zeros([0, 300, 300], dtype=np.int8)
img_mk_cat = np.zeros([0, 300, 300], dtype=np.int8)

for i in range(1, 2):
    filename = os.path.join('/Volumes/Macintosh HDD/CosasFer/test_data/Database_134_Angiograms', str(i+1) + '.pgm')
    print(filename)
    img = cv2.imread(filename, 0) / 255.0

    img_cat = np.concatenate((img_cat, img[np.newaxis, :, :]), axis=0)

start_time = time.time()
gmf_resp = gmf.gmfFilter(img_cat, 13, 9, np.array([2.0]), 12)
elapsed_time = time.time() - start_time
print('Computed in ', elapsed_time, ' seconds', gmf_resp.shape, np.max(gmf_resp), np.min(gmf_resp))

gmf_centered = 1.0 - (gmf_resp[:,0,:,:] - np.min(gmf_resp)) / (np.max(gmf_resp) - np.min(gmf_resp))

gmf_resp_2 = gmf.gmfFilter(gmf_centered, 13, 9, np.array([2.0]), 12)

gmf_resp_3 = gmf.gmfFilter(gmf_centered, 13, 9, np.array([0.05]), 12)

gmf_resp_4 = gmf.gmfFilter(gmf_centered, 13, 9, np.array([7.0]), 12)

print(gmf_resp_2.shape)
print(gmf_resp_3.shape)
print(gmf_resp_4.shape)

plt.subplot(2, 3, 2)
plt.imshow(gmf_resp[0,0])
plt.subplot(2, 3, 4)
plt.imshow(gmf_resp_2[0,0])
plt.subplot(2, 3, 5)
plt.imshow(gmf_resp_3[0,0])
plt.subplot(2, 3, 6)
plt.imshow(gmf_resp_4[0,0])
plt.show()