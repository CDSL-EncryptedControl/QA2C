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
    
from numpy.linalg import matrix_rank, eigvals
    
class arx:
    # sampling time
    ts = 0.01

    # observer controller state space model
    F = np.zeros((6,6), dtype=float)
    G = np.zeros((6,2), dtype=float)
    H = np.zeros((2,6), dtype=float)

    # gain for coordination transform
    L = np.zeros((6,2), dtype=float)

    # sequence of y and u like y0, y1, ...,y3 and u0, ...,u3
    # every y is column vector but save is row vector
    Ys = np.zeros((6,2), dtype=float)
    Us = np.zeros((6,2), dtype=float)

    # gain of y sequence and u sequence
    HG = np.zeros((12,2), dtype=float)
    HL = np.zeros((12,2), dtype=float)
    Href = np.zeros((12,6), dtype=float)

    # output
    u = np.zeros((2,1), dtype=float)

    def __init__(self, F, G, H, ts):
        # save original controller's state matrix
        self.ts = ts
        self.F = F
        self.G = G
        self.H = H
        
        O_matrix = ct.obsv(F, H)
        rank = matrix_rank(O_matrix)
        if rank != 6:
            is_obs = 0
            print(f"This controller is not observable")
        else:
            is_obs = 1
            
            # all desired_poles set to zero
            desired_poles = [0.0, 1e-5, -1e-5, 2e-6, -2e-6, 3e-7]

            # calc L, which is gain that F-LH can be nilpotent
            L = ct.place(self.F.T, self.H.T, desired_poles)

            # save gain L
            self.L = L.T
            
            # validation of gain L
            eigenvalues = eigvals(self.F - self.L @self.H)
            print(f"eigenvalue of F - LH:\n {eigenvalues}")

            # calc H(F-LH)^3G ... HG and H(F-LH)^3L ... HL
            self.HG[0:2,:] = (self.H @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ self.G)
            self.HG[2:4,:] = (self.H @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ self.G)
            self.HG[4:6,:] = (self.H @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ self.G)
            self.HG[6:8,:] = (self.H @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ self.G)
            self.HG[8:10,:] = (self.H @ (self.F - self.L @ self.H) @ self.G)    
            self.HG[10:12,:] = (self.H @ self.G)   

            self.HL[0:2,:] = (self.H @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ self.L)
            self.HL[2:4,:] = (self.H @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ self.L)
            self.HL[4:6,:] = (self.H @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ self.L)
            self.HL[6:8,:] = (self.H @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ self.L)
            self.HL[8:10,:] = (self.H @ (self.F - self.L @ self.H) @ self.L)    
            self.HL[10:12,:] = (self.H @ self.L)

            self.Href[0:2,:] = (self.H @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H))
            self.Href[2:4,:] = (self.H @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H))
            self.Href[4:6,:] = (self.H @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H))
            self.Href[6:8,:] = (self.H @ (self.F - self.L @ self.H) @ (self.F - self.L @ self.H))
            self.Href[8:10,:] = (self.H @ (self.F - self.L @ self.H))    
            self.Href[10:12,:] = (self.H)

    def mem_update(self, Yn, Un):
        # i = 0 is oldest value
        for i in range(5):
            self.Ys[i,:] = self.Ys[(i+1),:]
            self.Us[i,:] = self.Us[(i+1),:]

        # new value input in i = 5
        self.Ys[5,:] = Yn.reshape(1,2)
        self.Us[5,:] = Un.reshape(1,2)

    def get_output(self, ref):
        ref_concat = np.vstack((np.zeros((4,1), dtype=float), self.ts * ref))

        self.u[0, 0] = 0
        self.u[1, 0] = 0

        for i in range(6):
            self.u = self.u + self.HG[2*i:2*i+2,:] @ self.Ys[i,:].reshape(2,1) + self.HL[2*i:2*i+2,:] @ self.Us[i,:].reshape(2,1) + self.Href[2*i:2*i+2,:] @ ref_concat
        
        return self.u

class arx_q():
    ts = 0.01

    # quantized level s for matrix and r for signal
    s = 1
    r = 1

    # sequence of y and u like y0, y1, ...,y3 and u0, ...,u3
    # every y is column vector but save is row vector
    Ys_q = np.zeros((6,2), dtype=int)
    Us_q = np.zeros((6,2), dtype=int)

    # gain of y sequence and u sequence
    HG = np.zeros((12,2), dtype=float)
    HL = np.zeros((12,2), dtype=float)
    Href = np.zeros((12,6), dtype=float)
    HG_q = np.zeros((12,2), dtype=int)
    HL_q = np.zeros((12,2), dtype=int)
    Href_q = np.zeros((12,6), dtype=int)

    # output
    u_q = np.zeros((2,1), dtype=int)
    u = np.zeros((2,1), dtype=float)

    def __init__(self, HG, HL, Href, ts):
        self.HG = HG
        self.HL = HL
        self.Href = Href
        self.ts = ts

    def set_level(self, r, s):
        self.r = r
        self.s = s

    def quantize(self):
        self.HG_q = (self.s * self.HG).astype(int)
        self.HL_q = (self.s * self.HL).astype(int)
        self.Href_q = (self.s * self.Href).astype(int)

    def mem_update(self, Yn, Un):
        # i = 0 is oldest value
        for i in range(5):
            self.Ys_q[i,:] = self.Ys_q[(i+1),:]
            self.Us_q[i,:] = self.Us_q[(i+1),:]

        # new value input in i = 3
        self.Ys_q[5,:] = (self.r * Yn.reshape(1,2)).astype(int)
        self.Us_q[5,:] = (self.r * Un.reshape(1,2)).astype(int)

    def get_output(self, ref):
        ref_concat = np.vstack((np.zeros((4,1), dtype=float), self.ts * ref))
        ref_concat_q =  (self.r * ref_concat).astype(int)

        self.u_q[0,0] = 0
        self.u_q[1,0] = 0

        for i in range(6):
            self.u_q = self.u_q + self.HG_q[2*i:2*i+2,:] @ self.Ys_q[i,:].reshape(2,1) + self.HL_q[2*i:2*i+2,:] @ self.Us_q[i,:].reshape(2,1) + self.Href_q[2*i:2*i+2,:] @ ref_concat_q
        
        self.u[0, 0] = float(self.u_q[0, 0]) / self.r / self.s
        self.u[1, 0] = float(self.u_q[1, 0]) / self.r / self.s

        return self.u
