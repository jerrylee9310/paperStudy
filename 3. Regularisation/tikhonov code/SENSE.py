# -*- coding: utf-8 -*-

# =============================================================================
# # clear all; clc;
# from IPython import get_ipython
# get_ipython().magic('reset -sf') # '%reset -f' do the exactly same thing on the shell(non-interactive)
# get_ipython().magic('cls') # '%cls'
# =============================================================================

# check running time
import time

import scipy.io as spio
import numpy as np
import imgdomain
import kdomain


"""Tikhonov Regularisation"""
start_time = time.time()

# load k-space data
mat = spio.loadmat('Brain2D.mat', squeeze_me=True)
k_data = mat['DATA']
(ny,nx,nc) = np.shape(k_data)

# coil sensitivity
coil_sen = imgdomain.getcoilsen(kdomain.ktoimg(k_data))


# Parameters
r_fac = 4

# undersampling and corresponding image
under_k = kdomain.undersampling(k_data,r_fac)
under_img = kdomain.ktoimg(under_k)

# reconstruction
encd_func = np.zeros((nc,r_fac), dtype='complex')
recon_img = np.zeros((ny,nx), dtype='complex')
shifting_factor = ny//r_fac
vec_s = np.zeros((nc,1),dtype='float')

for x in range(nx):
    for y in range(shifting_factor):
        new_y = np.arange(y,ny,shifting_factor)
        
        for l in range(nc):
            encd_func[l,:] = (1/r_fac)*(coil_sen[new_y,x,l].reshape(1,r_fac))
        
        vec_s = under_img[y,x,:].reshape(nc,1)
        
        recon_img[new_y,x] = np.squeeze(np.matmul(np.linalg.pinv(encd_func),vec_s))
        


# show img
imgdomain.showimg(abs(recon_img))

print("--- %s seconds ---" % (time.time() - start_time))
    
    