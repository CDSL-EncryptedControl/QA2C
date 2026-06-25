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
    F = np.zeros((6,6), dtype=float)
    G = np.zeros((6,2), dtype=float)
    H = np.zeros((2,6), dtype=float)

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
        self.L = L.T

        # make a controller state space matrix
        F_horizon1 = np.hstack(((self.A - self.L @ self.C), np.zeros((4,2), dtype=float)))
        F_horizon2 = np.hstack(((np.zeros((2,4), dtype=float)), np.array([[1, 0], [0, 1]], dtype=float)))
        B_horizon = np.vstack((self.B, np.zeros((2,2), dtype=float)))
        self.F = np.vstack((F_horizon1, F_horizon2)) - B_horizon @ self.K
        self.G = np.vstack((self.L, np.array([[-self.ts, 0], [0, -self.ts]], dtype=float)))
        self.H = -self.K

    def state_update(self, y, u, ref):
        # update state on temp variable
        ref_concat = np.vstack((np.zeros((4,1), dtype=float), self.ts * ref))
        x_next = self.F @ self.x + self.G @ y + ref_concat
        
        # save state
        self.x = x_next
    
    def get_output(self):
        self.u = self.H @ self.x
        
        return self.u 
