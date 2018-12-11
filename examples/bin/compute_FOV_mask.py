#!/usr/bin/env python
# coding: utf-8

"""
CENTRO DE INVESTIGACION EN MATEMATICAS
DOCTORADO EN CIENCIAS DE LA COMPUTACION
FERNANDO CERVANTES SANCHEZ

FILE NAME: compute_FOV_mask.py

PURPOSE: Computes the FOV of the input image as a mask of 1's in the FOV and 0's in the rest.

FILE REFERENCES:
Name        I/O        Description
None        ----       ----------

ABNORMAL TERMINATION CONDITIONS, ERROR AND WARNING MESSAGES:

DEVELOPMENT HISTORY:
Date         Author        Change Id    Release    Description Of Change
20/Jun/2017  Fernando C.   0            1.0        Creation
"""

import numpy as np
import cv2
import matplotlib.pyplot as plt
import sys

def compute_FOV_mask(img, threshold_level = 0.1):
    height, width = img.shape
    mask = np.ones([height, width], dtype=np.uint8)
    threshold_value = (np.max(img) - np.min(img)) * threshold_level + np.min(img)
    mask[np.where(img < threshold_value)] = 0

    # erase the grouped pixels with size lesser than 1000 in the mask# Setup SimpleBlobDetector parameters.
    num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(mask, 8, cv2.CV_32S)
    
    for i in range(1, num_labels):
        if (stats[i][4] < 1000):
            labels[labels == i] = 0
    
    labels[labels > 0] = 255
    labels = labels.astype(np.uint8)

    # erode the resulting mask
    disk_mask = np.array([0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0], dtype=np.uint8).reshape([9, 9])
    mask = cv2.erode(labels, disk_mask)
    return mask
    

def mask_corners(mask):
    height, width, = mask.shape
    mask = 255 - mask
    
    # Check the labels of the corners pixels
    num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(mask, 8, cv2.CV_32S)
    
    # Check Upper Left corner:
    UL_label = labels[0, 0]
    
    # Check Upper Right corner:
    UR_label = labels[0, width - 1]
    
    # Check Lower Right corner:
    LR_label = labels[height - 1, width - 1]
    
    # Check Lower Left corner:
    LL_label = labels[height - 1, 0]
    
    masked_corners = np.zeros([height, width], dtype=np.uint8)
    
    if (UL_label > 0):
        masked_corners[labels == UL_label] = 255
    
    if (UR_label > 0):
        masked_corners[labels == UR_label] = 255
    
    if (LR_label > 0):
        masked_corners[labels == LR_label] = 255
    
    if (LL_label > 0):
        masked_corners[labels == LL_label] = 255
    
    return 255 - masked_corners
   

def fill_blackarea(img, mask, MAX_ITERS = 50):
    mask_temp = mask.copy()
    if np.max(img) > 1.0:
        img_temp = img / 255.0
    else:
        img_temp = img.copy()

    height, width, = img_temp.shape
    growing_kernel = np.array([0, 1,  0, 1, 1, 1, 0, 1, 0], dtype=np.uint8).reshape(3, 3)
    
    offset = np.arange(-10, 11)   
    
    growing_mask = mask_temp
    img_filled = img_temp
    
    n_iters = 1
    while(n_iters < MAX_ITERS):
        boundary = cv2.dilate(growing_mask, growing_kernel)
        boundary = np.logical_xor(boundary, growing_mask)

        if (np.sum(boundary) == 0):
            break
        
        y_bound, x_bound = np.where(boundary)
        
        for i in range(y_bound.shape[0]):
            if i == 0:
                print('iter: ' + str(n_iters) + ', border pixel: [' + str(x_bound[0]) + ', ' + str(y_bound[0]) + '], [' + str(x_bound[-1]) + ', ' + str(y_bound[-1]) + ']')
                
            x_rng = x_bound[i] + offset
            y_rng = y_bound[i] + offset

            x_rng = x_rng[x_rng >= 0]
            y_rng = y_rng[y_rng >= 0]
            
            x_rng = x_rng[x_rng < width]
            y_rng = y_rng[y_rng < height]
            
            active_mask = growing_mask[:, x_rng][y_rng, :]
            active_mask = active_mask.astype(np.bool)
            win_mean = img_filled[:, x_rng][y_rng, :]
            
            img_filled[y_bound[i], x_bound[i]] = np.sum(np.multiply(win_mean, active_mask)) / np.sum(active_mask)
        
        growing_mask = growing_mask.astype(np.bool)
        growing_mask = np.logical_or(growing_mask, boundary)
        growing_mask = growing_mask.astype(np.uint8)
                
        n_iters += 1
    
    return img_filled


if (__name__ == "__main__"):
    print("\tCompute FOV mask")
    # Read arguments:
    if (len(sys.argv) < 2):
        print("This program requires:")
        print('1) Path to an image')
        sys.exit()

    img = cv2.imread(sys.argv[1])
    img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    mask = compute_FOV_mask(img)
    mask = mask_corners(mask)
    img = fill_blackarea(img, mask)
        
    cv2.imshow("FOV mask", mask)
    cv2.waitKey(0)

    cv2.imshow("Filled image", img)
    cv2.waitKey(0)
