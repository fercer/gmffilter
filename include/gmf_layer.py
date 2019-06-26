import torch
import math

class ActiveArea(torch.autograd.Function):
    @staticmethod
    def forward(ctx, L, T, theta, delta, rows, cols):
        V = torch.linspace(-math.floor(rows/2.),math.floor(rows/2.)-1,rows)
        U = torch.linspace(-math.floor(cols/2.),math.floor(cols/2.)-1,cols)
        v, u = torch.meshgrid((V, U))
        
        cth = math.cos(theta)
        sth = math.sin(theta)
        
        # Compute the step functions
        left_step_fun_L = 1/(1+torch.exp(-(v*cth + u*sth + L/2.)*delta))
        right_step_fun_L = 1/(1+torch.exp(-(v*cth + u*sth - L/2.)*delta))
        
        left_step_fun_T = 1/(1+torch.exp(-(u*cth - v*sth + T/2.)*delta))
        right_step_fun_T = 1/(1+torch.exp(-(u*cth - v*sth - T/2.)*delta))
        
        h_v_L = left_step_fun_L - right_step_fun_L
        g_u_T = left_step_fun_T - right_step_fun_T
        
        # Define the box function from the active area is being computed
        box_fun_L_T = h_v_L * g_u_T;
        
        # Integrate along u and v axis
        box_fun_L_T_du = 0.5 * (box_fun_L_T[:,1:]+box_fun_L_T[:,:-1]).mv(U[1:]-U[:-1])
        box_fun_L_T_dudv = 0.5 * (box_fun_L_T_du[1:]+box_fun_L_T_du[:-1]).dot(V[1:]-V[:-1])
    
        ctx.save_for_backward(delta, left_step_fun_L, right_step_fun_L, left_step_fun_T, right_step_fun_T, h_v_L, g_u_T, U, V)
        
        return box_fun_L_T_dudv

    
    @staticmethod
    def backward(ctx, grad_output):
        delta, left_step_fun_L, right_step_fun_L, left_step_fun_T, right_step_fun_T, h_v_L, g_u_T, U, V = ctx.saved_tensors
        dleft_step_fun_L = left_step_fun_L * (1-left_step_fun_L) * delta/2.
        dright_step_fun_L = right_step_fun_L * (1-right_step_fun_L) * delta/2.
        dh_v_L = dleft_step_fun_L + dright_step_fun_L
        
        dleft_step_fun_T = left_step_fun_T * (1-left_step_fun_T) * delta/2.
        dright_step_fun_T = right_step_fun_T * (1-right_step_fun_T) * delta/2.
        dg_u_T = dleft_step_fun_T + dright_step_fun_T
        
        dbox_fun_dL_T = dh_v_L * g_u_T;
        dbox_fun_L_dT = h_v_L * dg_u_T;
        
        # Integrate along u and v axis the partial derivative w.r.t. L
        dbox_fun_dL_T_du = 0.5 * (dbox_fun_dL_T[:,1:]+dbox_fun_dL_T[:,:-1]).mv(U[1:]-U[:-1])
        dbox_fun_dL_T_dudv = 0.5 * (dbox_fun_dL_T_du[1:]+dbox_fun_dL_T_du[:-1]).dot(V[1:]-V[:-1])
        
        # Integrate along u and v axis the partial derivative w.r.t. T
        dbox_fun_L_dT_du = 0.5 * (dbox_fun_L_dT[:,1:]+dbox_fun_L_dT[:,:-1]).mv(U[1:]-U[:-1])
        dbox_fun_L_dT_dudv = 0.5 * (dbox_fun_L_dT_du[1:]+dbox_fun_L_dT_du[:-1]).dot(V[1:]-V[:-1])
        
        return dbox_fun_dL_T_dudv.view(1), dbox_fun_L_dT_dudv.view(1), None, None, None, None


class ProfileArea(torch.autograd.Function):
    @staticmethod
    def forward(ctx, sigma, L, T, theta, delta, rows, cols):
        V = torch.linspace(-math.floor(rows/2.),math.floor(rows/2.)-1,rows)
        U = torch.linspace(-math.floor(cols/2.),math.floor(cols/2.)-1,cols)
        v, u = torch.meshgrid((V, U))
        
        cth = math.cos(theta)
        sth = math.sin(theta)
        
        # Compute the profile function
        f_s = torch.exp(-(u*cth - v*sth)**2/(2*sigma**2))
        
        # Compute the step functions
        left_step_fun_L = 1/(1+torch.exp(-(v*cth + u*sth + L/2.)*delta))
        right_step_fun_L = 1/(1+torch.exp(-(v*cth + u*sth - L/2.)*delta))
        
        left_step_fun_T = 1/(1+torch.exp(-(u*cth - v*sth + T/2.)*delta))
        right_step_fun_T = 1/(1+torch.exp(-(u*cth - v*sth - T/2.)*delta))
        
        h_v_L = left_step_fun_L - right_step_fun_L
        g_u_T = left_step_fun_T - right_step_fun_T
        
        # Define the box function from the active area is being computed
        box_fun_s_L_T = f_s * h_v_L * g_u_T;
        
        # Integrate along u and v axis
        box_fun_s_L_T_du = 0.5 * (box_fun_s_L_T[:,1:]+box_fun_s_L_T[:,:-1]).mv(U[1:]-U[:-1])
        box_fun_s_L_T_dudv = 0.5 * (box_fun_s_L_T_du[1:]+box_fun_s_L_T_du[:-1]).dot(V[1:]-V[:-1])
    
        ctx.save_for_backward(sigma, theta, delta, f_s, left_step_fun_L, right_step_fun_L, left_step_fun_T, right_step_fun_T, h_v_L, g_u_T, U, V, u, v)
        
        return box_fun_s_L_T_dudv

    
    @staticmethod
    def backward(ctx, grad_output):
        sigma, theta, delta, f_s, left_step_fun_L, right_step_fun_L, left_step_fun_T, right_step_fun_T, h_v_L, g_u_T, U, V, u, v = ctx.saved_tensors
        df_s = (u*math.cos(theta) - v*math.sin(theta))**2/sigma**3 * f_s
        
        dleft_step_fun_L = left_step_fun_L * (1-left_step_fun_L) * delta/2.
        dright_step_fun_L = right_step_fun_L * (1-right_step_fun_L) * delta/2.
        dh_v_L = dleft_step_fun_L + dright_step_fun_L
        
        dleft_step_fun_T = left_step_fun_T * (1-left_step_fun_T) * delta/2.
        dright_step_fun_T = right_step_fun_T * (1-right_step_fun_T) * delta/2.
        dg_u_T = dleft_step_fun_T + dright_step_fun_T
        
        dbox_fun_ds_L_T = df_s * h_v_L * g_u_T;
        dbox_fun_s_dL_T = f_s * dh_v_L * g_u_T;
        dbox_fun_s_L_dT = f_s * h_v_L * dg_u_T;
        
        # Integrate along u and v axis the partial derivative w.r.t. sigma
        dbox_fun_ds_L_T_du = 0.5 * (dbox_fun_ds_L_T[:,1:]+dbox_fun_ds_L_T[:,:-1]).mv(U[1:]-U[:-1])
        dbox_fun_ds_L_T_dudv = 0.5 * (dbox_fun_ds_L_T_du[1:]+dbox_fun_ds_L_T_du[:-1]).dot(V[1:]-V[:-1])
                
        # Integrate along u and v axis the partial derivative w.r.t. L
        dbox_fun_s_dL_T_du = 0.5 * (dbox_fun_s_dL_T[:,1:]+dbox_fun_s_dL_T[:,:-1]).mv(U[1:]-U[:-1])
        dbox_fun_s_dL_T_dudv = 0.5 * (dbox_fun_s_dL_T_du[1:]+dbox_fun_s_dL_T_du[:-1]).dot(V[1:]-V[:-1])
        
        # Integrate along u and v axis the partial derivative w.r.t. T
        dbox_fun_s_L_dT_du = 0.5 * (dbox_fun_s_L_dT[:,1:]+dbox_fun_s_L_dT[:,:-1]).mv(U[1:]-U[:-1])
        dbox_fun_s_L_dT_dudv = 0.5 * (dbox_fun_s_L_dT_du[1:]+dbox_fun_s_L_dT_du[:-1]).dot(V[1:]-V[:-1])
        
        return dbox_fun_ds_L_T_dudv.view(1), dbox_fun_s_dL_T_dudv.view(1), dbox_fun_s_L_dT_dudv.view(1), None, None, None, None


if __name__ == '__main__':
    print('Testing active area and profile area for generation of a smooth GMF filter')
    sigma = torch.autograd.Variable(torch.Tensor([2.0]), requires_grad=True)
    L = torch.autograd.Variable(torch.Tensor([9]), requires_grad=True)
    T = torch.autograd.Variable(torch.Tensor([13]), requires_grad=True)
    theta = torch.autograd.Variable(torch.Tensor([10./180.*math.pi]), requires_grad=False)
    delta = torch.autograd.Variable(torch.Tensor([1.]), requires_grad=False)
    
    my_active_area = ActiveArea.apply
    operation = my_active_area(L, T, theta, delta, 512, 512)
    operation.backward()
    print('dActive_area/dL', L.grad)
    print('dActive_area/dT', T.grad)

    L.grad.data.zero_()
    T.grad.data.zero_()

    my_profile_area = ProfileArea.apply
    operation = my_profile_area(sigma, L, T, theta, delta, 512, 512)
    operation.backward()
    print('dProfile_area/dsigma', sigma.grad)
    print('dProfile_area/dL', L.grad)
    print('dProfile_area/dT', T.grad)
