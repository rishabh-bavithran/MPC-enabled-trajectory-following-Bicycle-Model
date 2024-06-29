import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
import numpy as np

m=1500
Iz=3000
Caf=19000
Car=33000
lf=2
lr=3
Ts=0.02

Q=np.matrix('35 0;0 1') #1-4 step weights
S=np.matrix('35 0;0 1') #final step weights
R=np.matrix('1') # input weight

outputs=2 
hz = 20
x_dot=20
lane_width=7
nr_lanes=5
r=4
f=0.01
time_length = 5

constants={'m':m, 'Iz':Iz, 'Caf':Caf, 'Car':Car, 'lf':lf, 'lr':lr, 'Ts':Ts, 'Q':Q, 'S':S, 'R':R, 'outputs':outputs, 'hz':hz, 'x_dot':x_dot, 'r':r, 'f':f, 'time_length':time_length, 'lane_width':lane_width}

def trajectory_generator(t):

    x_dot=constants['x_dot']
    x=np.linspace(0,x_dot*t[-1],num=len(t))
    y = 8 * np.sin((np.pi/6) * t) - 8 * np.cos((np.pi/6) * t) 

    dx=x[1:len(x)]-x[0:len(x)-1]
    dy=y[1:len(y)]-y[0:len(y)-1]
    psi=np.zeros(len(x))
    psiInt=psi
    psi[0]=np.arctan2(dy[0],dx[0])
    psi[1:len(psi)]=np.arctan2(dy[0:len(dy)],dx[0:len(dx)])

    dpsi=psi[1:len(psi)]-psi[0:len(psi)-1]
    psiInt[0]=psi[0]
    for i in range(1,len(psiInt)):
        if dpsi[i-1]<-np.pi:
            psiInt[i]=psiInt[i-1]+(dpsi[i-1]+2*np.pi)
        elif dpsi[i-1]>np.pi:
            psiInt[i]=psiInt[i-1]+(dpsi[i-1]-2*np.pi)
        else:
            psiInt[i]=psiInt[i-1]+dpsi[i-1]
    return psiInt,x,y

def state_space():
    m=constants['m']
    Iz=constants['Iz']
    Caf=constants['Caf']
    Car=constants['Car']
    lf=constants['lf']
    lr=constants['lr']
    Ts=constants['Ts']
    x_dot=constants['x_dot']

    A1=-(2*Caf+2*Car)/(m*x_dot)
    A2=-x_dot-(2*Caf*lf-2*Car*lr)/(m*x_dot)
    A3=-(2*lf*Caf-2*lr*Car)/(Iz*x_dot)
    A4=-(2*lf**2*Caf+2*lr**2*Car)/(Iz*x_dot)

    A=np.array([[A1, 0, A2, 0],[0, 0, 1, 0],[A3, 0, A4, 0],[1, x_dot, 0, 0]])
    B=np.array([[2*Caf/m],[0],[2*lf*Caf/Iz],[0]])
    C=np.array([[0, 1, 0, 0],[0, 0, 0, 1]])
    D=0

    Ad=np.identity(np.size(A,1))+Ts*A
    Bd=Ts*B
    Cd=C
    Dd=D

    return Ad, Bd, Cd, Dd

def mpc_simplification(Ad, Bd, Cd, Dd, hz):

    A_aug=np.concatenate((Ad,Bd),axis=1)
    temp1=np.zeros((np.size(Bd,1),np.size(Ad,1)))
    temp2=np.identity(np.size(Bd,1))
    temp=np.concatenate((temp1,temp2),axis=1)

    A_aug=np.concatenate((A_aug,temp),axis=0)
    B_aug=np.concatenate((Bd,np.identity(np.size(Bd,1))),axis=0)
    C_aug=np.concatenate((Cd,np.zeros((np.size(Cd,0),np.size(Bd,1)))),axis=1)
    D_aug=Dd

    Q=constants['Q']
    S=constants['S']
    R=constants['R']

    CQC=np.matmul(np.transpose(C_aug),Q)
    CQC=np.matmul(CQC,C_aug)

    CSC=np.matmul(np.transpose(C_aug),S)
    CSC=np.matmul(CSC,C_aug)

    QC=np.matmul(Q,C_aug)
    SC=np.matmul(S,C_aug)


    Qdb=np.zeros((np.size(CQC,0)*hz,np.size(CQC,1)*hz))
    Tdb=np.zeros((np.size(QC,0)*hz,np.size(QC,1)*hz))
    Rdb=np.zeros((np.size(R,0)*hz,np.size(R,1)*hz))
    Cdb=np.zeros((np.size(B_aug,0)*hz,np.size(B_aug,1)*hz))
    Adc=np.zeros((np.size(A_aug,0)*hz,np.size(A_aug,1)))

    for i in range(0,hz):
        if i == hz-1:
            Qdb[np.size(CSC,0)*i:np.size(CSC,0)*i+CSC.shape[0],np.size(CSC,1)*i:np.size(CSC,1)*i+CSC.shape[1]]=CSC
            Tdb[np.size(SC,0)*i:np.size(SC,0)*i+SC.shape[0],np.size(SC,1)*i:np.size(SC,1)*i+SC.shape[1]]=SC
        else:
            Qdb[np.size(CQC,0)*i:np.size(CQC,0)*i+CQC.shape[0],np.size(CQC,1)*i:np.size(CQC,1)*i+CQC.shape[1]]=CQC
            Tdb[np.size(QC,0)*i:np.size(QC,0)*i+QC.shape[0],np.size(QC,1)*i:np.size(QC,1)*i+QC.shape[1]]=QC

        Rdb[np.size(R,0)*i:np.size(R,0)*i+R.shape[0],np.size(R,1)*i:np.size(R,1)*i+R.shape[1]]=R

        for j in range(0,hz):
            if j<=i:
                Cdb[np.size(B_aug,0)*i:np.size(B_aug,0)*i+B_aug.shape[0],np.size(B_aug,1)*j:np.size(B_aug,1)*j+B_aug.shape[1]]=np.matmul(np.linalg.matrix_power(A_aug,((i+1)-(j+1))),B_aug)

        Adc[np.size(A_aug,0)*i:np.size(A_aug,0)*i+A_aug.shape[0],0:0+A_aug.shape[1]]=np.linalg.matrix_power(A_aug,i+1)

    Hdb=np.matmul(np.transpose(Cdb),Qdb)
    Hdb=np.matmul(Hdb,Cdb)+Rdb

    temp=np.matmul(np.transpose(Adc),Qdb)
    temp=np.matmul(temp,Cdb)

    temp2=np.matmul(-Tdb,Cdb)
    Fdbt=np.concatenate((temp,temp2),axis=0)

    return Hdb,Fdbt,Cdb,Adc

def open_loop_new_states(states,U1):

    m=constants['m']
    Iz=constants['Iz']
    Caf=constants['Caf']
    Car=constants['Car']
    lf=constants['lf']
    lr=constants['lr']
    Ts=constants['Ts']
    x_dot=constants['x_dot']

    current_states=states
    new_states=current_states
    y_dot=current_states[0]
    psi=current_states[1]
    psi_dot=current_states[2]
    Y=current_states[3]

    sub_loop=30
    for i in range(0,sub_loop):
        y_dot_dot=-(2*Caf+2*Car)/(m*x_dot)*y_dot+(-x_dot-(2*Caf*lf-2*Car*lr)/(m*x_dot))*psi_dot+2*Caf/m*U1
        psi_dot=psi_dot
        psi_dot_dot=-(2*lf*Caf-2*lr*Car)/(Iz*x_dot)*y_dot-(2*lf**2*Caf+2*lr**2*Car)/(Iz*x_dot)*psi_dot+2*lf*Caf/Iz*U1
        Y_dot=np.sin(psi)*x_dot+np.cos(psi)*y_dot
        y_dot=y_dot+y_dot_dot*Ts/sub_loop
        psi=psi+psi_dot*Ts/sub_loop
        psi_dot=psi_dot+psi_dot_dot*Ts/sub_loop
        Y=Y+Y_dot*Ts/sub_loop

    new_states[0]=y_dot
    new_states[1]=psi
    new_states[2]=psi_dot
    new_states[3]=Y

    return new_states

constants=constants
Ts=constants['Ts']
outputs=constants['outputs']
hz = constants['hz'] 
x_dot=constants['x_dot']
time_length=constants['time_length']

t=np.arange(0,time_length+Ts,Ts)
r=constants['r']
f=constants['f']
psi_ref,X_ref,Y_ref=trajectory_generator(t)
sim_length=len(t) 
refSignals=np.zeros(len(X_ref)*outputs)
k=0
for i in range(0,len(refSignals),outputs):
    refSignals[i]=psi_ref[k]
    refSignals[i+1]=Y_ref[k]
    k=k+1

y_dot=0.
psi=0.
psi_dot=0.
Y=Y_ref[0]-3.
states=np.array([y_dot,psi,psi_dot,Y])
statesTotal=np.zeros((len(t),len(states))) 
statesTotal[0][0:len(states)]=states
psi_opt_total=np.zeros((len(t),hz))
Y_opt_total=np.zeros((len(t),hz))
U1=0 
UTotal=np.zeros(len(t))
UTotal[0]=U1
C_psi_opt=np.zeros((hz,(len(states)+np.size(U1))*hz))
for i in range(1,hz+1):
    C_psi_opt[i-1][i+4*(i-1)]=1

C_Y_opt=np.zeros((hz,(len(states)+np.size(U1))*hz))
for i in range(3,hz+3):
    C_Y_opt[i-3][i+4*(i-3)]=1

Ad,Bd,Cd,Dd=state_space()
Hdb,Fdbt,Cdb,Adc=mpc_simplification(Ad,Bd,Cd,Dd,hz)
k=0
for i in range(0,sim_length-1):
    x_aug_t=np.transpose([np.concatenate((states,[U1]),axis=0)])
    k=k+outputs
    if k+outputs*hz<=len(refSignals):
        r=refSignals[k:k+outputs*hz]
    else:
        r=refSignals[k:len(refSignals)]
        hz=hz-1

    if hz<constants['hz']: 
        Hdb,Fdbt,Cdb,Adc=mpc_simplification(Ad,Bd,Cd,Dd,hz)

    ft=np.matmul(np.concatenate((np.transpose(x_aug_t)[0][0:len(x_aug_t)],r),axis=0),Fdbt)
    du=-np.matmul(np.linalg.inv(Hdb),np.transpose([ft]))
    x_aug_opt=np.matmul(Cdb,du)+np.matmul(Adc,x_aug_t)
    psi_opt=np.matmul(C_psi_opt[0:hz,0:(len(states)+np.size(U1))*hz],x_aug_opt)
    Y_opt=np.matmul(C_Y_opt[0:hz,0:(len(states)+np.size(U1))*hz],x_aug_opt)
    psi_opt=np.transpose((psi_opt))[0]
    psi_opt_total[i+1][0:hz]=psi_opt
    Y_opt=np.transpose((Y_opt))[0]
    Y_opt_total[i+1][0:hz]=Y_opt
    U1=U1+du[0][0]
    if U1 < -np.pi/6:
        U1=-np.pi/6
    elif U1 > np.pi/6:
        U1=np.pi/6
    else:
        U1=U1
    UTotal[i+1]=U1
    states=open_loop_new_states(states,U1)
    statesTotal[i+1][0:len(states)]=states

frame_amount=int(time_length/Ts)
lf=constants['lf']
lr=constants['lr']
def update_plot(num):
    hz = constants['hz'] 
    car_1.set_data([X_ref[num]-lr*np.cos(statesTotal[num,1]),X_ref[num]+lf*np.cos(statesTotal[num,1])],
        [statesTotal[num,3]-lr*np.sin(statesTotal[num,1]),statesTotal[num,3]+lf*np.sin(statesTotal[num,1])])
    if num+hz>len(t):
        hz=len(t)-num
    if num!=0:
        car_predicted.set_data(X_ref[num:num+hz],Y_opt_total[num][0:hz])
    car_determined.set_data(X_ref[0:num],statesTotal[0:num,3])
    return car_1, car_determined, car_predicted,

fig_x=15
fig_y=5
fig=plt.figure(figsize=(fig_x,fig_y),dpi=120,facecolor=(0.8,0.8,0.8))
n=1
m=1
gs=gridspec.GridSpec(n,m)

ax0=fig.add_subplot(gs[:,:],facecolor=(0.2,0.2,0.2))
ref_trajectory=ax0.plot(X_ref,Y_ref,'g', linestyle='--', linewidth=1, label = "Reference Trajectory")
lane_1,=ax0.plot([X_ref[0],X_ref[frame_amount]],[0,0],color = (0.7,0.7,0.7), linestyle='--',linewidth=6, label= "Lane Divider")
car_1,=ax0.plot([],[],'w',linewidth=3, label= "Bicycle")
car_predicted,=ax0.plot([],[],color = 'orange',linewidth=1, label="Horizon Period Prediction")
car_determined,=ax0.plot([],[],'-r',linewidth=1, label= "Traversed Trajectory")
plt.xlim(X_ref[0],X_ref[frame_amount])
plt.ylim(-X_ref[frame_amount]/(n*(fig_x/fig_y)*2),X_ref[frame_amount]/(n*(fig_x/fig_y)*2))
plt.ylabel('Lanes',fontsize=15)
plt.xlabel('distance travelled [m]',fontsize=15)
ax0.legend(loc='lower right', fontsize='large')
car_ani=animation.FuncAnimation(fig, update_plot,
    frames=frame_amount,interval=20,repeat=True,blit=True)
plt.gcf().canvas.set_window_title('MPC_Bicycle_Model')
plt.show()
