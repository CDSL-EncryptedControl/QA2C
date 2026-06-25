import numpy as np
import control as ct

class intc:
    # sampling period
    ts = 0.1

    # parameters
    Jp = 0.0211;  Jy = 0.0221
    Dp = 0.0053;  Dy = 0.0062
    Mg = 0.0153
    Kpp =  0.0011;  Kyy =  0.0047
    Kpy =  0.0021;  Kyp = -0.0027

    # system(state) matrix(discrete model)
    A = np.zeros((4,4), dtype=float)
    B = np.zeros((4,2), dtype=float)
    C = np.zeros((2,4), dtype=float)
    D = np.zeros((2,2), dtype=float)

    # output gain K and observer gain L
    K = np.zeros((2,6), dtype=float)
    L = np.zeros((4,2), dtype=float)

    # observer controller state space model
    F = np.zeros((4,4), dtype=float)
    G = np.zeros((4,2), dtype=float)
    H = np.zeros((2,4), dtype=float)
    J = np.zeros((2,2), dtype=float) 

    # state and output
    x = np.zeros((6,1), dtype=float)
    u = np.zeros((2,1), dtype=float)

    def __init__(self, ts):
        # write your continous time linear model
        A = np.array([[0, 0, 1, 0],
                      [0, 0, 0, 1],
                      [-self.Mg/self.Jp, 0, -self.Dp/self.Jp, 0],
                      [0, 0, 0, -self.Dy/self.Jy]], dtype=float)
        B = np.array([[0, 0],
                      [0, 0],
                      [self.Kpp/self.Jp, self.Kpy/self.Jp],
                      [self.Kyp/self.Jy, self.Kyy/self.Jy]], dtype=float)
        C = np.array([[1, 0, 0, 0],
                      [0, 1, 0, 0]], dtype=float)
        D = np.array([[0, 0],
                      [0, 0]], dtype=float)
        
        A_horizon = np.hstack((A, np.zeros((4,2), dtype=float)))
        C_horizon = np.hstack((-C, np.zeros((2,2), dtype=float)))
        integral_A = np.vstack((A_horizon, C_horizon))
        integral_B = np.vstack((B, np.zeros((2,2), dtype=float)))
        integral_C = np.hstack((C, np.zeros((2,2), dtype=float)))
        sys_c = ct.ss(A, B, C, D)
        sys_d = sys_c.sample(ts, method='zoh')

        sys_c_int = ct.ss(integral_A, integral_B, integral_C, D)
        sys_d_int = sys_c_int.sample(ts, method='zoh')

        # save discretization value
        self.A = sys_d.A
        self.B = sys_d.B
        self.C = sys_d.C
        self.D = D
        self.ts = ts

        # for gain K dlqr parameters setting
        Q_k = np.array([[1500, 0, 0, 0, 0, 0],
                        [0, 1500, 0, 0, 0, 0],
                        [0, 0, 10, 0, 0, 0],
                        [0, 0, 0, 10, 0, 0],
                        [0, 0, 0, 0, 200, 0],
                        [0, 0, 0, 0, 0, 200]], dtype=float)
        R_k = np.array([[0.1, 0],
                        [0, 0.1]], dtype=float)
        K, Sk, Ek = ct.dlqr(sys_d_int.A, sys_d_int.B, Q_k, R_k)

        # for gain L dlqr parameters setting
        Q_l = np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 100, 0],
                        [0, 0, 0, 100]], dtype=float)
        R_l = np.array([[0.1, 0],
                        [0, 0.1]], dtype=float)
        L, Sl, El = ct.dlqr(self.A.T, self.C.T, Q_l, R_l)

        # # for gain L pole-place parameters setting
        # L = ct.place(self.A.T, self.C.T, [0.65, 0.66, 0.67, 0.68])

        # save gain K and N_bar
        self.K = K
        self.L[0:4,:] = L.T

        # make a controller state space matrix
        self.F = self.A
        self.G = self.L
        self.H = -self.K

    def state_update(self, y, u, ref):
        # update state on temp variable
        x_obs = self.F @ self.x[0:4,:] + self.B @ u + self.G @ (y - self.C @ self.x[0:4,:])
        x_int = self.x[4:6,:] + self.ts * (ref - y)
        x_next = np.vstack((x_obs, x_int))
        
        # save state
        self.x = x_next
    
    def get_output(self):
        self.u = self.H @ self.x
        
        return self.u
        
