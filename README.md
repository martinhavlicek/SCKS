# SCKS toolbox

Implementaion of Square-root Cubature Kalman Filter& Square-root Rauch-Tang-Striebel Smoother for modeling time sequences as described in [1,2]. Tested with Matlab R2012b.

This algorithm performs joint estimation of the states, input and parameters of the model that is described as a stochastic continuous-discrete state-space in terms of nonlinear blind deconvolution. The state equations must have the form of ordinary differential equations, where the discretization is performed through local-linearization scheme [3]. Additionally, the parameter noise covariance is estimated online via stochastic Robbins-Monro approximation method [4], and the measurement noise covariance is estimated online as well by using combination of varitional Bayesian (VB) approach with nonlinear filter/smoother [5].

### Demos:
* *demo_CKS_convolution.m*: Example of linear convolution system inversion.
* *demo_CKS_doublewell.m*: Example of (nonlinear) double-well system inversion.
* *demo_CKS_lorenz.m*: Example of lorenz system inversion.
* *demo_CKS_hdm.m*: Example of hemodynamic model inversion.
* *demo_CKS_hdm_HR.m*:  Example of hemodynamic model inversion in high-temporal-resolution mode.

### References:
1. Havlicek et al. (2011) Modeling neuronal responses in fMRI using cubature Kalman filter. Neuroimage, [doi](https://doi.org/10.1016/j.neuroimage.2011.03.005).
2. Arasaratnam, I., Haykin, S. (2009) Cubature Kalman Filters. IEEE Transactions on Automatic Control 54, 1254-1269.
3. Jimenez, J.C. (2002) A simple algebraic expression to evaluate the local linearization schemes for stochastic differential equations, Applied Mathematics Letters 15, 775-780.
4. Van der Merwe, R., 2004. Sigma-point Kalman filters for probabilistic inference in dynamic state-space models. Ph.D.thesis, Oregon Graduate Institute of Science and Technology.
5. Sarkka, S., Hartikainen, J. (2013) Extension of VB-AKF to Estimation of Full Covariance and Non-Linear Systems. arXiv:1302.0681v1.


