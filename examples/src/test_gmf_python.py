import cv2
import numpy as np
import matplotlib.pyplot as plt
import gmf

import compute_FOV_mask as fov
import computeAUC as roc

import os
import sys
import time

print('Testing GMF filter')

filename = os.path.join(r'E:/test_data/Angios_134', '17.png')    
print(filename)

img = np.zeros([512, 512])
img[106:-106,106:-106] = cv2.imread(filename, 0) / 255.0

mask = fov.compute_FOV_mask(img)
img = fov.fill_blackarea(img, mask)

filename = os.path.join(r'E:/test_data/Angios_134', '17_gt.png')
img_gt = cv2.imread(filename, 0)
img_gt = img_gt / np.max(img_gt)

print('Img GT:', img_gt.dtype)

fft_img = np.fft.fft2(np.fft.fftshift(img))

resp_dog = np.zeros([12, 512, 512])
resp_gmf = np.zeros([12, 512, 512])

U, V = np.meshgrid(np.arange(-256,256), np.arange(-256,256))

test_sigma = 2.65
test_L = 13
test_T = 6.0 * test_sigma
alpha = 1.0 + 1.0/(6.0 * test_sigma * test_L)

for k, test_theta in enumerate(np.linspace(0.0, np.pi - np.pi/12.0, 12)):
    ctheta = np.cos(test_theta)
    stheta = np.sin(test_theta)
    u = U*ctheta - V*stheta
    v = V*ctheta + U*stheta

    step_upper = 1.0 / (1.0 + np.exp(-1e2*(v+test_L/2.0))) 
    step_lower = 1.0 / (1.0 + np.exp(-1e2*(v-test_L/2.0))) 
    step_V = step_upper - step_lower

    step_upper = 1.0 / (1.0 + np.exp(-1e2*(u+test_T/2.0)))
    step_lower = 1.0 / (1.0 + np.exp(-1e2*(u-test_T/2.0)))
    step_H = step_upper - step_lower

    gaussian = np.exp(-0.5 * u**2 / test_sigma**2)

    area_A = np.sum(step_H * step_V)
    area_P = np.sum(step_H * step_V * gaussian)

    GMF = ((1.0 - gaussian) / (area_A - area_P) - 1.0/area_A) * step_H * step_V

    alpha_gaussian = alpha * np.exp(-0.5*alpha**2*u**2/test_sigma**2)
    DOG = gaussian - alpha_gaussian
    gaussian_A = np.sum(gaussian[U.shape[0]//2,:])
    gaussian_B = np.sum(alpha_gaussian[U.shape[0]//2,:])
    print('Area A: {}/{}'.format(gaussian_A, np.sqrt(2.0*np.pi)*test_sigma), 'area B: {}/{}'.format(gaussian_B, alpha*np.sqrt(2.0*np.pi)*test_sigma/alpha))

    DOG = DOG * step_V

    print('DOG sum', np.sum(DOG))
    print('GMF sum', np.sum(GMF), 'profile area:', area_P, 'area:', area_A)

    fft_dog = np.fft.fft2(DOG)
    fft_gmf = np.fft.fft2(GMF)

    fft_resp_dog = fft_img * fft_dog
    fft_resp_gmf = fft_img * fft_gmf

    resp_dog[k,:,:] = np.fft.ifft2(fft_resp_dog).real
    resp_gmf[k,:,:] = np.fft.ifft2(fft_resp_gmf).real

resp_dog = np.max(resp_dog, axis=0)[106:-106,106:-106]
resp_gmf = np.max(resp_gmf, axis=0)[106:-106,106:-106]
resp_gmf_o = gmf.gmfFilter(img, test_L, np.array([test_sigma]), 12, 1)[0,0][106:-106,106:-106]

print('Diff:', np.sqrt(np.sum((resp_gmf - resp_dog)**2)))

plt.subplot(1, 4, 1)
plt.imshow(resp_dog)
plt.subplot(1, 4, 2)
plt.imshow(resp_gmf)
plt.subplot(1, 4, 3)
plt.imshow(resp_gmf_o)
plt.subplot(1, 4, 4)
plt.imshow(img_gt)
plt.show()

Az_dog = roc.aucROC(resp_dog, img_gt)
Az_gmf = roc.aucROC(resp_gmf, img_gt)
Az_gmf_o = roc.aucROC(resp_gmf_o, img_gt)
print('DOG roc auc: {}, GMF roc auc: {}, Org roc auc: {}'.format(Az_dog, Az_gmf, Az_gmf_o))

"""
start_time = time.time()
gmf_resp = gmf.gmfFilter(img_cat, 9, np.array([2.0]), 12, 1)
elapsed_time = time.time() - start_time

print(gmf_resp.shape, gmf_resp.max())
print('Computed in ', elapsed_time, ' seconds', gmf_resp.shape, np.max(gmf_resp), np.min(gmf_resp))
"""
"""
gmf_centered = 1.0 - (gmf_resp[:,0,:,:] - np.min(gmf_resp)) / (np.max(gmf_resp) - np.min(gmf_resp))

gmf_resp_2 = gmf.gmfFilter(gmf_centered, 9, np.array([2.0]), 12)

gmf_resp_3 = gmf.gmfFilter(gmf_centered, 9, np.array([0.05]), 12)

gmf_resp_4 = gmf.gmfFilter(gmf_centered, 9, np.array([7.0]), 12)

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
"""
"""
gmf_resp_K = gmf.gmfFilter(img_cat, 9, np.array([2.0]), 12, 0)
print('Max resp:', gmf_resp_K.max())
print(np.sqrt(np.sum(gmf_resp_K[0].max(axis=0) - gmf_resp[0, 0])**2))

plt.imshow(gmf_resp[0, 0])
plt.show()
"""