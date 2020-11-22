import numpy as np

def undersampling(k_data,reduction_factor):
    """ Undersampling along y-direction """
    
    
    under_k = np.zeros_like(k_data, dtype='complex')
    
        
    (ny,nx,nc) = np.shape(k_data)
    
    for slice_idx in range(nc):
        data_slice = np.squeeze(k_data[:,:,slice_idx])
        
        # shift zero frequency to centre
        kshift = np.fft.fftshift(data_slice)
        
        # initialize a subsampled array with complex numbers
        subshift = np.zeros_like(kshift)
        
        for i in range(0,ny,reduction_factor):
            subshift[i] = kshift[i]
        
        under_k[:,:,slice_idx] = np.fft.fftshift(subshift)
        
    return under_k


def ktoimg(k_data):
    (ny,nx,nc) = np.shape(k_data)
    
    under_img = np.zeros_like(k_data, dtype='complex')

    for slice_idx in range(nc):
        data_slice = np.squeeze(k_data[:,:,slice_idx])
        
        # shift zero frequency to centre
        kshift = np.fft.fftshift(data_slice)
        
        under_img[:,:,slice_idx] = np.fft.ifftshift((np.fft.ifft2(kshift)).astype(complex))
        
    return under_img