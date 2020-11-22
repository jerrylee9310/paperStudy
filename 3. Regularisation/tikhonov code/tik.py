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
import matplotlib.pyplot as plt
import imgdomain
import kdomain


"""Tikhonov Regularisation"""
start_time = time.time()

# load k-space data
mat = spio.loadmat('Brain2D.mat', squeeze_me=True)
k_data = mat['DATA']
(ny,nx,nc) = np.shape(k_data)

# coil sensitivity
full_img = kdomain.ktoimg(k_data)
coil_sen = imgdomain.getcoilsen(full_img)

# inter-channel coil noise covariance
phi = imgdomain.intercoil_noise_covariance(full_img,25)
phi_inv = np.linalg.inv(phi)  
w,v = np.linalg.eig(phi_inv)  # phi_inv = v @ np.diag(w) @ np.matrix.getH(v) // diagonalisation

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


# g_factor map initialisation
g_map = np.zeros((ny,nx), dtype='complex')
brutal_iter = 500

x_zero = np.zeros((r_fac,1))
err_model = np.zeros((brutal_iter,1), dtype='float')
err_prior = np.zeros((brutal_iter,1), dtype='float')

# regularisation parameter // 10**(-7) ~ 10**(-8)
for k in range(brutal_iter):
    lmbd = (1+k)*10**(-9)
    
    for x in range(nx):
        for y in range(shifting_factor):
            new_y = np.arange(y,ny,shifting_factor)
            
            for l in range(nc):
                encd_func[l,:] = (1/r_fac)*(coil_sen[new_y,x,l].reshape(1,r_fac))
            
            vec_s = under_img[y,x,:].reshape(nc,1)
            
            # noise whitening // A tilde, y tilde
            A_tilde = np.diag(w**(-1/2)) @ np.matrix.getH(v) @ encd_func # encoding matrix
            y_tilde = np.diag(w**(-1/2)) @ np.matrix.getH(v) @ vec_s # aliased image
            
            U,S,Vh = np.linalg.svd(A_tilde)
            V = np.matrix.getH(Vh)
            
            
            
            
            # tikhonov regularisation
            large_gamma = np.concatenate( ( np.diag(S/(S**2 + lmbd**2)) , np.zeros((r_fac,nc-r_fac), dtype='float') ) , axis = 1)
            large_phi = np.diag(lmbd**2/(S**2 + lmbd**2))
            
            x_lambda = V @ large_gamma @ np.matrix.getH(U) @ y_tilde + V @ large_phi @ Vh @ x_zero
            
            # error
            err_model[k,0] = err_model[k,0] + np.sum((abs(y_tilde - A_tilde @ x_lambda))**2, axis = 0)
            err_prior[k,0] = err_prior[k,0] + np.sum((abs(x_lambda - x_zero))**2, axis = 0)
            
    
    print(k)
            
find_lmbd = err_model + err_prior
find_lmbd_scaling = 10**(14)*err_model + err_prior
min_lmbd = np.ndarray.tolist(find_lmbd_scaling)
min_lmbd_idx = min_lmbd.index(min(min_lmbd))

fig1, ax1 = plt.subplots()
ax1.plot(find_lmbd,'--', label='model_err + prior_err')
ax1.plot(find_lmbd_scaling, label = '10^14*model_err + prior_err')
ax1.plot(min_lmbd_idx,min_lmbd[min_lmbd_idx],'ro', markersize = 8, label = 'min')
legend = ax1.legend(loc='upper center', shadow=True, fontsize='x-large')
plt.show()

fig2, ax2 = plt.subplots()
ax2.plot(err_prior, err_model)
ax2.plot(err_prior[min_lmbd_idx],err_model[min_lmbd_idx],'ro')
ax2.set_xlabel('Prior error')
ax2.set_ylabel('Model error')
plt.show()

    



lmbd = (1+min_lmbd_idx)*10**(-9)
    
for x in range(nx):
    for y in range(shifting_factor):
        new_y = np.arange(y,ny,shifting_factor)
        
        for l in range(nc):
            encd_func[l,:] = (1/r_fac)*(coil_sen[new_y,x,l].reshape(1,r_fac))
        
        vec_s = under_img[y,x,:].reshape(nc,1)
        
        # noise whitening // A tilde, y tilde
        A_tilde = np.diag(w**(-1/2)) @ np.matrix.getH(v) @ encd_func # encoding matrix
        y_tilde = np.diag(w**(-1/2)) @ np.matrix.getH(v) @ vec_s # aliased image
        
        U,S,Vh = np.linalg.svd(A_tilde)
        V = np.matrix.getH(Vh)
        
        # tikhonov regularisation
        large_gamma = np.concatenate( ( np.diag(S/(S**2 + lmbd**2)) , np.zeros((r_fac,nc-r_fac), dtype='float') ) , axis = 1)
        large_phi = np.diag(lmbd**2/(S**2 + lmbd**2))
        
        x_lambda = V @ large_gamma @ np.matrix.getH(U) @ y_tilde + V @ large_phi @ Vh @ x_zero
        
        # g_facotr
        unfolded_var = np.linalg.inv(V @ large_gamma @ np.matrix.getH(large_gamma) @ Vh)
        
        full_var = np.zeros((r_fac,1), dtype = 'complex')
        for rf in range(r_fac):
            U_ful, S_ful, Vh_ful = np.linalg.svd(A_tilde[:,rf].reshape((nc,1)))
            full_var[rf,0] = np.matrix.getH(Vh_ful) @ np.diag(S_ful**2) @ Vh_ful
         
        for g_idx in range(r_fac):
            g_map[new_y[g_idx],x] = np.sqrt(unfolded_var[g_idx][g_idx] * full_var[g_idx,0])
                
        recon_img[new_y,x] = np.squeeze(np.matmul(np.linalg.pinv(encd_func),vec_s))

        

# show img
plt.subplot(121); plt.imshow(abs(recon_img),'gray'); plt.title('recon img')
plt.colorbar()
plt.subplot(122); plt.imshow(abs(g_map),'gray'); plt.title('g-factor')
plt.colorbar()

print("--- %s seconds ---" % (time.time() - start_time))
    
    