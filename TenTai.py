#This class is to create the function of tidal-evolution

import numpy as np
import scipy.integrate

class TenTai:
    G=6.674e-11

    def __init__ (self,Mp,Ms,
                      Rp,Rs,
                      Dtp,#Dts,
                      k2p,k2s,
                      Qp,#Qs
                      ADt,
                      AQ
                      ):
        self.Mp=Mp
        self.Ms=Ms
        self.Rp=Rp
        self.Rs=Rs
        self.Dtp=Dtp
#        self.Dts=Dts
        self.k2p=k2p
        self.k2s=k2s
        self.Qp=Qp
#        self.Qs=Qs
        self.ADt=ADt
        self.AQ=AQ
    def mu(self):
        return self.Mp*self.Ms/(self.Mp+self.Ms)
    def Cp(self):
        return 0.328*self.Mp*self.Rp**2
        #return 0.3308*self.Mp*self.Rp**2
    def Cs(self):
        return 0.4*self.Ms*self.Rs**2
        #return 0.3929*self.Mp*self.Rp**2
    def Dts(self):
        return self.Dtp*self.ADt*(self.k2p/self.k2s)*(self.Ms/self.Mp)**2*(self.Rp/self.Rs)**5
    def Qs(self):
        return self.Qp/self.AQ*(self.k2s/self.k2p)*(self.Mp/self.Ms)**2*(self.Rs/self.Rp)**5
    def Dt_model(self,t,x):
        Eqs= np.zeros((4))
        Eqs[0]= (-3*TenTai.G/(self.Cp()*x[2]**6)*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-1)*((1+15/2*x[3]**2)*x[0]-(1+27/2*x[3]**2))
                +3/2*x[0]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*((1+27/2*x[3]**2)*(x[0]+self.ADt*x[1])-(1+23*x[3]**2)*(1+self.ADt)))*31536000
        Eqs[1]= (-3*TenTai.G/(self.Cs()*x[2]**6)*self.k2s*self.Dts()*self.Mp**2*(self.Rs**5/self.Rp**6)*((1+15/2*x[3]**2)*x[1]-(1+27/2*x[3]**2))
                +3/2*x[1]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*((1+27/2*x[3]**2)*(x[0]+self.ADt*x[1])-(1+23*x[3]**2)*(1+self.ADt)))*31536000
        Eqs[2]= (x[2]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*((1+27/2*x[3]**2)*(x[0]+self.ADt*x[1])-(1+23*x[3]**2)*(1+self.ADt)))*31536000
        Eqs[3]= (x[3]*(27*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*((11/18*(x[0]+self.ADt*x[1])-(1+self.ADt))))*31536000   
        
        return Eqs
    
    def Dt_model_exactly(self,t,x):
        Eqs= np.zeros((4))
        Eqs[0]= (-3*TenTai.G/(self.Cp()*x[2]**6)*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-1)*(((1+3*x[3]**2+3/8*x[3]**4)/(1-x[3]**2)**(9/2))*x[0]-((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)/(1-x[3]**2)**6))
                +3/2*x[0]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*(((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)/(1-x[3]**2)**6)*(x[0]+self.ADt*x[1])-((1+31/2*x[3]**2+255/8*x[3]**4+185/16*x[3]**6+25/64*x[3]**8)/(1-x[3]**2)**7.5)*(1+self.ADt)))*31536000
        Eqs[1]= (-3*TenTai.G/(self.Cs()*x[2]**6)*self.k2s*self.Dts()*self.Mp**2*(self.Rs**5/self.Rp**6)*(((1+3*x[3]**2+3/8*x[3]**4)/(1-x[3]**2)**(9/2))*x[1]-((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)/(1-x[3]**2)**6))
                +3/2*x[1]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*(((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)/(1-x[3]**2)**6)*(x[0]+self.ADt*x[1])-((1+31/2*x[3]**2+255/8*x[3]**4+185/16*x[3]**6+25/64*x[3]**8)/(1-x[3]**2)**7.5)*(1+self.ADt)))*31536000
        Eqs[2]= (x[2]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*(((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)/(1-x[3]**2)**6)*(x[0]+self.ADt*x[1])-((1+31/2*x[3]**2+255/8*x[3]**4+185/16*x[3]**6+25/64*x[3]**8)/(1-x[3]**2)**7.5)*(1+self.ADt)))*31536000
        Eqs[3]= (x[3]*(27*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*((((1+3/2*x[3]**2+1/8*x[3]**4)/(1-x[3]**2)**5)*11/18*(x[0]+self.ADt*x[1])-((1+15/4*x[3]**2+15/8*x[3]**4+5/64*x[3]**6)/(1-x[3]**2)**6.5)*(1+self.ADt))))*31536000   
        
        return Eqs
    
    
    def Q_model(self,t,x,Dp,Ep,Fp,Ds,Es,Fs,Sp,Ss):
        Eqs= np.zeros((4)) 
        Eqs[0]= (-3*(TenTai.G)/((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))/(2*self.Cp()*x[2]**6)*self.k2p/self.Qp*self.Ms**2*(self.Rp**-1)*(Sp+x[3]**2*Dp)
                +3/2*x[0]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(Sp+self.AQ*Ss+x[3]**2*(Ep+self.AQ*Es))))*31536000
        Eqs[1]=(-3*(TenTai.G)/((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))/(2*self.Cs()*x[2]**6)*self.k2s/self.Qs()*self.Mp**2*(self.Rs**5/self.Rp**6)*(Ss+x[3]**2*Ds)
                +3/2*x[1]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(Sp+self.AQ*Ss+x[3]**2*(Ep+self.AQ*Es))))*31536000
        Eqs[2]= x[2]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(Sp+self.AQ*Ss+x[3]**2*(Ep+self.AQ*Es)))*31536000
        Eqs[3]= x[3]*((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(Fp+self.AQ*Fs))*31536000
          
        return Eqs  
    
    def Q_model_exactly(self,t,x,Dp,Ep,Fp,Ds,Es,Fs,Sp,Ss,Gp,Hp,Jp,Gs,Hs,Js):
        Eqs= np.zeros((4)) 
        Eqs[0]= (-3*(TenTai.G)/((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))/(2*self.Cp()*x[2]**6)*self.k2p/self.Qp*self.Ms**2*(self.Rp**-1)*(Sp+x[3]**2*Dp+x[3]**4*Gp)
                +3/2*x[0]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(Sp+self.AQ*Ss+x[3]**2*(Ep+self.AQ*Es)+x[3]**4*(Hp+self.AQ*Hs))))*31536000
        Eqs[1]=(-3*(TenTai.G)/((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))/(2*self.Cs()*x[2]**6)*self.k2s/self.Qs()*self.Mp**2*(self.Rs**5/self.Rp**6)*(Ss+x[3]**2*Ds+x[3]**4*Gs)
                +3/2*x[1]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(Sp+self.AQ*Ss+x[3]**2*(Ep+self.AQ*Es)+x[3]**4*(Hp+self.AQ*Hs))))*31536000
        Eqs[2]= x[2]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(Sp+self.AQ*Ss+x[3]**2*(Ep+self.AQ*Es)+x[3]**4*(Hp+self.AQ*Hs)))*31536000
        Eqs[3]= x[3]*((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(Fp+self.AQ*Fs+x[3]**2*(Jp+self.AQ*Js)))*31536000
          
        return Eqs  
    
    def Q_model_exactly_continuous(self,t,x):
              
        b=1000000
        
        Eqs= np.zeros((4)) 
        Eqs[0]= (-3*(TenTai.G)/((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))/(2*self.Cp()*x[2]**6)*self.k2p/self.Qp*self.Ms**2*(self.Rp**-1)*(np.tanh(b*(2*x[0]-2))+x[3]**2*(1/4*np.tanh(b*(1*x[0]-2))-5*np.tanh(b*(2*x[0]-2))+49/4*np.tanh(b*(2*x[0]-3)))+x[3]**4*(-1/16*np.tanh(b*(2*x[0]-1))+63/8*np.tanh(b*(2*x[0]-2))-861/16*np.tanh(b*(2*x[0]-3))+289/4*np.tanh(b*(2*x[0]-4))))
                +3/2*x[0]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(np.tanh(b*(2*x[0]-2))+self.AQ*np.tanh(b*(2*x[1]-2))+x[3]**2*((3/4*np.tanh(b*(-1))+1/8*np.tanh(b*(2*x[0]-1))-5*np.tanh(b*(2*x[0]-2))+147/8*np.tanh(b*(2*x[0]-3)))+self.AQ*(3/4*np.tanh(b*(-1))+1/8*np.tanh(b*(2*x[1]-1))-5*np.tanh(b*(2*x[1]-2))+147/8*np.tanh(b*(2*x[1]-3))))+x[3]**4*((32/27*np.tanh(b*(-1))+81/16*np.tanh(b*(-2))+115/64*np.tanh(b*(2*x[0]-1))-6*np.tanh(b*(2*x[0]-2))-3579/64*np.tanh(b*(2*x[0]-3))-527/4*np.tanh(b*(2*x[0]-4)))+self.AQ*(32/27*np.tanh(b*(-1))+81/16*np.tanh(b*(-2))+115/64*np.tanh(b*(2*x[1]-1))-6*np.tanh(b*(2*x[1]-2))-3579/64*np.tanh(b*(2*x[1]-3))-527/4*np.tanh(b*(2*x[1]-4)))))))*31536000
        Eqs[1]=(-3*(TenTai.G)/((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))/(2*self.Cs()*x[2]**6)*self.k2s/self.Qs()*self.Mp**2*(self.Rs**5/self.Rp**6)*(np.tanh(b*(2*x[1]-2))+x[3]**2*(1/4*np.tanh(b*(2*x[1]-1))-5*np.tanh(b*(2*x[1]-2))+49/4*np.tanh(b*(2*x[1]-3)))+x[3]**4*(-1/16*np.tanh(b*(2*x[1]-1))+63/8*np.tanh(b*(2*x[1]-2))-861/16*np.tanh(b*(2*x[1]-3))+289/4*np.tanh(b*(2*x[1]-4))))
                +3/2*x[0]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(np.tanh(b*(2*x[0]-2))+self.AQ*np.tanh(b*(2*x[1]-2))+x[3]**2*((3/4*np.tanh(b*(-1))+1/8*np.tanh(b*(2*x[0]-1))-5*np.tanh(b*(2*x[0]-2))+147/8*np.tanh(b*(2*x[0]-3)))+self.AQ*(3/4*np.tanh(b*(-1))+1/8*np.tanh(b*(2*x[1]-1))-5*np.tanh(b*(2*x[1]-2))+147/8*np.tanh(b*(2*x[1]-3))))+x[3]**4*((32/27*np.tanh(b*(-1))+81/16*np.tanh(b*(-2))+115/64*np.tanh(b*(2*x[0]-1))-6*np.tanh(b*(2*x[0]-2))-3579/64*np.tanh(b*(2*x[0]-3))-527/4*np.tanh(b*(2*x[0]-4)))+self.AQ*(32/27*np.tanh(b*(-1))+81/16*np.tanh(b*(-2))+115/64*np.tanh(b*(2*x[1]-1))-6*np.tanh(b*(2*x[1]-2))-3579/64*np.tanh(b*(2*x[1]-3))-527/4*np.tanh(b*(2*x[1]-4)))))))*31536000
        Eqs[2]= x[2]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(np.tanh(b*(2*x[0]-2))+self.AQ*np.tanh(b*(2*x[1]-2))+x[3]**2*((3/4*np.tanh(b*(-1))+1/8*np.tanh(b*(2*x[0]-1))-5*np.tanh(b*(2*x[0]-2))+147/8*np.tanh(b*(2*x[0]-3)))+self.AQ*(3/4*np.tanh(b*(-1))+1/8*np.tanh(b*(2*x[1]-1))-5*np.tanh(b*(2*x[1]-2))+147/8*np.tanh(b*(2*x[1]-3))))+x[3]**4*((32/27*np.tanh(b*(-1))+81/16*np.tanh(b*(-2))+115/64*np.tanh(b*(2*x[0]-1))-6*np.tanh(b*(2*x[0]-2))-3579/64*np.tanh(b*(2*x[0]-3))-527/4*np.tanh(b*(2*x[0]-4)))+self.AQ*(32/27*np.tanh(b*(-1))+81/16*np.tanh(b*(-2))+115/64*np.tanh(b*(2*x[1]-1))-6*np.tanh(b*(2*x[1]-2))-3579/64*np.tanh(b*(2*x[1]-3))-527/4*np.tanh(b*(2*x[1]-4))))))*31536000
        Eqs[3]= 3/4*x[3]*((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*((3/2*np.tanh(b*(-1))-1/4*np.tanh(b*(2*x[0]-1))-1*np.tanh(b*(2*x[0]-2))+49/4*np.tanh(b*(2*x[0]-3)))+self.AQ*(3/2*np.tanh(b*(-1))-1/4*np.tanh(b*(2*x[1]-1))-1*np.tanh(b*(2*x[1]-2))+49/4*np.tanh(b*(2*x[1]-3)))+x[3]**2*((3/16*np.tanh(b*(-1))+81/8*np.tanh(b*(-2))+119/32*np.tanh(b*(2*x[0]-1))-45/2*np.tanh(b*(2*x[0]-2))-919/32*np.tanh(b*(2*x[0]-3))+119*np.tanh(b*(2*x[0]-4)))+self.AQ*(3/16*np.tanh(b*(-1))+81/8*np.tanh(b*(-2))+119/32*np.tanh(b*(2*x[1]-1))-45/2*np.tanh(b*(2*x[1]-2))-919/32*np.tanh(b*(2*x[1]-3))+119*np.tanh(b*(2*x[1]-4))))))*31536000
        
        return Eqs 
    
    def Q_model_exactly_continuous4(self,t,x):
              
        b=1000000
        
        Eqs= np.zeros((4)) 
        Eqs[0]= (-3*(TenTai.G)/((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))/(2*self.Cp()*x[2]**6)*self.k2p/self.Qp*self.Ms**2*(self.Rp**-1)*(np.tanh(b*(2*x[0]-2)))
                +3/2*x[0]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*((np.tanh(b*(2*x[0]-2)))+self.AQ*(np.tanh(b*(2*x[1]-2))))))*31536000
        Eqs[1]=(-3*(TenTai.G)/((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))/(2*self.Cs()*x[2]**6)*self.k2s/self.Qs()*self.Mp**2*(self.Rs**5/self.Rp**6)*(np.tanh(b*(2*x[1]-2)))
                +3/2*x[0]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*((np.tanh(b*(2*x[0]-2)))+self.AQ*(np.tanh(b*(2*x[1]-2))))))*31536000
        Eqs[2]= x[2]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*((np.tanh(b*(2*x[0]-2)))+self.AQ*(np.tanh(b*(2*x[1]-2)))))*31536000
        Eqs[3]= 0
         
        return Eqs
        
    def Q_model_exactly_continuous2(self,t,x):
              
        b=100000000
        
        Eqs= np.zeros((4)) 
        Eqs[0]= (-3*(TenTai.G)/((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))/(2*self.Cp()*x[2]**6)*self.k2p/self.Qp*self.Ms**2*(self.Rp**-1)*(np.tanh(4*b*x[1]/(1-(b*x[1]**2)))+x[3]**2*(1/4*np.tanh(b*(1*x[0]-2))-5*np.tanh(b*(2*x[0]-2))+49/4*np.tanh(b*(2*x[0]-3)))+x[3]**4*(-1/16*np.tanh(b*(2*x[0]-1))+63/8*np.tanh(b*(2*x[0]-2))-861/16*np.tanh(b*(2*x[0]-3))+289/4*np.tanh(b*(2*x[0]-4))))
                +3/2*x[0]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(np.tanh(4*b*x[1]/(1-(b*x[1]**2)))+self.AQ*np.tanh(b*(2*x[1]-2))+x[3]**2*((3/4*np.tanh(b*(-1))+1/8*np.tanh(b*(2*x[0]-1))-5*np.tanh(b*(2*x[0]-2))+147/8*np.tanh(b*(2*x[0]-3)))+self.AQ*(3/4*np.tanh(b*(-1))+1/8*np.tanh(b*(2*x[1]-1))-5*np.tanh(b*(2*x[1]-2))+147/8*np.tanh(b*(2*x[1]-3))))+x[3]**4*((32/27*np.tanh(b*(-1))+81/16*np.tanh(b*(-2))+115/64*np.tanh(b*(2*x[0]-1))-6*np.tanh(b*(2*x[0]-2))-3579/64*np.tanh(b*(2*x[0]-3))-527/4*np.tanh(b*(2*x[0]-4)))+self.AQ*(32/27*np.tanh(b*(-1))+81/16*np.tanh(b*(-2))+115/64*np.tanh(b*(2*x[1]-1))-6*np.tanh(b*(2*x[1]-2))-3579/64*np.tanh(b*(2*x[1]-3))-527/4*np.tanh(b*(2*x[1]-4)))))))*31536000
        Eqs[1]=(-3*(TenTai.G)/((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))/(2*self.Cs()*x[2]**6)*self.k2s/self.Qs()*self.Mp**2*(self.Rs**5/self.Rp**6)*(np.tanh(b*(2*x[1]-2))+x[3]**2*(1/4*np.tanh(b*(2*x[1]-1))-5*np.tanh(b*(2*x[1]-2))+49/4*np.tanh(b*(2*x[1]-3)))+x[3]**4*(-1/16*np.tanh(b*(2*x[1]-1))+63/8*np.tanh(b*(2*x[1]-2))-861/16*np.tanh(b*(2*x[1]-3))+289/4*np.tanh(b*(2*x[1]-4))))
                +3/2*x[0]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(np.tanh(4*b*x[1]/(1-(b*x[1]**2)))+self.AQ*np.tanh(b*(2*x[1]-2))+x[3]**2*((3/4*np.tanh(b*(-1))+1/8*np.tanh(b*(2*x[0]-1))-5*np.tanh(b*(2*x[0]-2))+147/8*np.tanh(b*(2*x[0]-3)))+self.AQ*(3/4*np.tanh(b*(-1))+1/8*np.tanh(b*(2*x[1]-1))-5*np.tanh(b*(2*x[1]-2))+147/8*np.tanh(b*(2*x[1]-3))))+x[3]**4*((32/27*np.tanh(b*(-1))+81/16*np.tanh(b*(-2))+115/64*np.tanh(b*(2*x[0]-1))-6*np.tanh(b*(2*x[0]-2))-3579/64*np.tanh(b*(2*x[0]-3))-527/4*np.tanh(b*(2*x[0]-4)))+self.AQ*(32/27*np.tanh(b*(-1))+81/16*np.tanh(b*(-2))+115/64*np.tanh(b*(2*x[1]-1))-6*np.tanh(b*(2*x[1]-2))-3579/64*np.tanh(b*(2*x[1]-3))-527/4*np.tanh(b*(2*x[1]-4)))))))*31536000
        Eqs[2]= x[2]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(np.tanh(4*b*x[1]/(1-(b*x[1]**2)))+self.AQ*np.tanh(b*(2*x[1]-2))+x[3]**2*((3/4*np.tanh(b*(-1))+1/8*np.tanh(b*(2*x[0]-1))-5*np.tanh(b*(2*x[0]-2))+147/8*np.tanh(b*(2*x[0]-3)))+self.AQ*(3/4*np.tanh(b*(-1))+1/8*np.tanh(b*(2*x[1]-1))-5*np.tanh(b*(2*x[1]-2))+147/8*np.tanh(b*(2*x[1]-3))))+x[3]**4*((32/27*np.tanh(b*(-1))+81/16*np.tanh(b*(-2))+115/64*np.tanh(b*(2*x[0]-1))-6*np.tanh(b*(2*x[0]-2))-3579/64*np.tanh(b*(2*x[0]-3))-527/4*np.tanh(b*(2*x[0]-4)))+self.AQ*(32/27*np.tanh(b*(-1))+81/16*np.tanh(b*(-2))+115/64*np.tanh(b*(2*x[1]-1))-6*np.tanh(b*(2*x[1]-2))-3579/64*np.tanh(b*(2*x[1]-3))-527/4*np.tanh(b*(2*x[1]-4))))))*31536000
        Eqs[3]= 3/4*x[3]*((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*((3/2*np.tanh(b*(-1))-1/4*np.tanh(b*(2*x[0]-1))-1*np.tanh(b*(2*x[0]-2))+49/4*np.tanh(b*(2*x[0]-3)))+self.AQ*(3/2*np.tanh(b*(-1))-1/4*np.tanh(b*(2*x[1]-1))-1*np.tanh(b*(2*x[1]-2))+49/4*np.tanh(b*(2*x[1]-3)))+x[3]**2*((3/16*np.tanh(b*(-1))+81/8*np.tanh(b*(-2))+119/32*np.tanh(b*(2*x[0]-1))-45/2*np.tanh(b*(2*x[0]-2))-919/32*np.tanh(b*(2*x[0]-3))+119*np.tanh(b*(2*x[0]-4)))+self.AQ*(3/16*np.tanh(b*(-1))+81/8*np.tanh(b*(-2))+119/32*np.tanh(b*(2*x[1]-1))-45/2*np.tanh(b*(2*x[1]-2))-919/32*np.tanh(b*(2*x[1]-3))+119*np.tanh(b*(2*x[1]-4))))))*31536000
        
        return Eqs 
    def Q_model_exactly_continuous3(self,t,x):
              
        b=100000
        
        Eqs= np.zeros((4)) 
        Eqs[0]= (-3*(TenTai.G)/((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))/(2*self.Cp()*x[2]**6)*self.k2p/self.Qp*self.Ms**2*(self.Rp**-1)*(np.tanh(b*(2*x[0]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+x[3]**2*(1/4*np.tanh(b*(1*x[0]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-5*np.tanh(b*(2*x[0]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+49/4*np.tanh(b*(2*x[0]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))))+x[3]**4*(-1/16*np.tanh(b*(2*x[0]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+63/8*np.tanh(b*(2*x[0]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-861/16*np.tanh(b*(2*x[0]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+289/4*np.tanh(b*(2*x[0]-4*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))))
                +3/2*x[0]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(np.tanh(b*(2*x[0]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+self.AQ*np.tanh(b*(2*x[1]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+x[3]**2*((3/4*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+1/8*np.tanh(b*(2*x[0]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-5*np.tanh(b*(2*x[0]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+147/8*np.tanh(b*(2*x[0]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))))+self.AQ*(3/4*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+1/8*np.tanh(b*(2*x[1]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-5*np.tanh(b*(2*x[1]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+147/8*np.tanh(b*(2*x[1]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))))+x[3]**4*((32/27*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+81/16*np.tanh(b*(-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+115/64*np.tanh(b*(2*x[0]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-6*np.tanh(b*(2*x[0]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-3579/64*np.tanh(b*(2*x[0]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-527/4*np.tanh(b*(2*x[0]-4*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))))+self.AQ*(32/27*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+81/16*np.tanh(b*(-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+115/64*np.tanh(b*(2*x[1]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-6*np.tanh(b*(2*x[1]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-3579/64*np.tanh(b*(2*x[1]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-527/4*np.tanh(b*(2*x[1]-4*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))))))))*31536000
        Eqs[1]=(-3*(TenTai.G)/((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))/(2*self.Cs()*x[2]**6)*self.k2s/self.Qs()*self.Mp**2*(self.Rs**5/self.Rp**6)*(np.tanh(b*(2*x[1]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+x[3]**2*(1/4*np.tanh(b*(2*x[1]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-5*np.tanh(b*(2*x[1]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+49/4*np.tanh(b*(2*x[1]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))))+x[3]**4*(-1/16*np.tanh(b*(2*x[1]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+63/8*np.tanh(b*(2*x[1]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-861/16*np.tanh(b*(2*x[1]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+289/4*np.tanh(b*(2*x[1]-4*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))))
                +3/2*x[0]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(np.tanh(b*(2*x[0]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+self.AQ*np.tanh(b*(2*x[1]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+x[3]**2*((3/4*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+1/8*np.tanh(b*(2*x[0]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-5*np.tanh(b*(2*x[0]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+147/8*np.tanh(b*(2*x[0]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))))+self.AQ*(3/4*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+1/8*np.tanh(b*(2*x[1]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-5*np.tanh(b*(2*x[1]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+147/8*np.tanh(b*(2*x[1]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))))+x[3]**4*((32/27*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+81/16*np.tanh(b*(-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+115/64*np.tanh(b*(2*x[0]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-6*np.tanh(b*(2*x[0]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-3579/64*np.tanh(b*(2*x[0]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-527/4*np.tanh(b*(2*x[0]-4*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))))+self.AQ*(32/27*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+81/16*np.tanh(b*(-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+115/64*np.tanh(b*(2*x[1]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-6*np.tanh(b*(2*x[1]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-3579/64*np.tanh(b*(2*x[1]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-527/4*np.tanh(b*(2*x[1]-4*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))))))))*31536000
        Eqs[2]= x[2]*(3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*(np.tanh(b*(2*x[0]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+self.AQ*np.tanh(b*(2*x[1]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+x[3]**2*((3/4*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+1/8*np.tanh(b*(2*x[0]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-5*np.tanh(b*(2*x[0]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+147/8*np.tanh(b*(2*x[0]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))))+self.AQ*(3/4*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+1/8*np.tanh(b*(2*x[1]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-5*np.tanh(b*(2*x[1]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+147/8*np.tanh(b*(2*x[1]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))))+x[3]**4*((32/27*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+81/16*np.tanh(b*(-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+115/64*np.tanh(b*(2*x[0]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-6*np.tanh(b*(2*x[0]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-3579/64*np.tanh(b*(2*x[0]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-527/4*np.tanh(b*(2*x[0]-4*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))))+self.AQ*(32/27*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+81/16*np.tanh(b*(-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+115/64*np.tanh(b*(2*x[1]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-6*np.tanh(b*(2*x[1]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-3579/64*np.tanh(b*(2*x[1]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-527/4*np.tanh(b*(2*x[1]-4*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))))))*31536000
        Eqs[3]= 3/4*x[3]*((TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)*self.k2p/self.Qp*(self.Ms/self.Mp)/x[2]**5*((3/2*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-1/4*np.tanh(b*(2*x[0]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-1*np.tanh(b*(2*x[0]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+49/4*np.tanh(b*(2*x[0]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))))+self.AQ*(3/2*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-1/4*np.tanh(b*(2*x[1]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-1*np.tanh(b*(2*x[1]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+49/4*np.tanh(b*(2*x[1]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))))+x[3]**2*((3/16*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+81/8*np.tanh(b*(-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+119/32*np.tanh(b*(2*x[0]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-45/2*np.tanh(b*(2*x[0]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-919/32*np.tanh(b*(2*x[0]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+119*np.tanh(b*(2*x[0]-4*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5))))+self.AQ*(3/16*np.tanh(b*(-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+81/8*np.tanh(b*(-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+119/32*np.tanh(b*(2*x[1]-1*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-45/2*np.tanh(b*(2*x[1]-2*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))-919/32*np.tanh(b*(2*x[1]-3*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))+119*np.tanh(b*(2*x[1]-4*(TenTai.G*(self.Mp+self.Ms)/(x[2]*self.Rp)**3)**(0.5)))))))*31536000
          
        return Eqs 
    
    def Dt_model_change_w(self,t,x):
        Eqs= np.zeros((4))
        Eqs[0]= (-3*TenTai.G/(self.Cp()*x[2]**6)*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-1)*((1+15/2*x[3]**2)*(x[0]+0.999)-(1+27/2*x[3]**2))
                +3/2*x[0]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*((1+27/2*x[3]**2)*((x[0]+0.999)+self.ADt*(x[1]+0.999))-(1+23*x[3]**2)*(1+self.ADt)))*31536000
        Eqs[1]= (-3*TenTai.G/(self.Cs()*x[2]**6)*self.k2s*self.Dts()*self.Mp**2*(self.Rs**5/self.Rp**6)*((1+15/2*x[3]**2)*(x[1]+0.999)-(1+27/2*x[3]**2))
                +3/2*x[1]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*((1+27/2*x[3]**2)*((x[0]+0.999)+self.ADt*(x[1]+0.999))-(1+23*x[3]**2)*(1+self.ADt)))*31536000
        Eqs[2]= (x[2]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*((1+27/2*x[3]**2)*((x[0]+0.999)+self.ADt*(x[1]+0.999))-(1+23*x[3]**2)*(1+self.ADt)))*31536000
        Eqs[3]= (x[3]*(27*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*((11/18*((x[0]+0.999)+self.ADt*(x[1]+0.999))-(1+self.ADt))))*31536000   
        
        return Eqs


    def Dt_model_nots(self,t,x):
        Eqs= np.zeros((4))
        Eqs[0]= (-3*TenTai.G/(self.Cp()*x[2]**6)*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-1)*((1+15/2*x[3]**2)*(x[0]+1)-(1+27/2*x[3]**2))
                +3/2*x[0]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*((1+27/2*x[3]**2)*(x[0]+1)-(1+23*x[3]**2)*1))*31536000
        Eqs[1]= 0
        Eqs[2]= (x[2]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*((1+27/2*x[3]**2)*(x[0]+1)-(1+23*x[3]**2)*1))*31536000
        Eqs[3]= (x[3]*(27*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*((11/18*(x[0]+1)-1)))*31536000   
        
        return Eqs
     
    
#    def i_Dt_model(self,t,x):
#        Eqs= np.zeros((5))
#        Eqs[0]= (-3*TenTai.G/(self.Cp()*x[2]**6)*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-1)*(((1+3*x[3]**2+3/8*x[3]**4)/(1-x[3]**2)**(9/2))*x[0]*np.cos(x[4])-((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)/(1-x[3]**2)**6))
#                +3/2*x[0]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*(((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)/(1-x[3]**2)**6)*(x[0]*np.cos(x[4])+self.ADt*x[1])-((1+31/2*x[3]**2+255/8*x[3]**4+185/16*x[3]**6+25/64*x[3]**8)/(1-x[3]**2)**7.5)*(1+self.ADt)))*31536000
#        Eqs[1]= (-3*TenTai.G/(self.Cs()*x[2]**6)*self.k2s*self.Dts()*self.Mp**2*(self.Rs**5/self.Rp**6)*(((1+3*x[3]**2+3/8*x[3]**4)/(1-x[3]**2)**(9/2))*x[1]*np.cos(x[4])-((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)/(1-x[3]**2)**6))
#                +3/2*x[0]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*(((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)/(1-x[3]**2)**6)*(x[0]*np.cos(x[4])+self.ADt*x[1])-((1+31/2*x[3]**2+255/8*x[3]**4+185/16*x[3]**6+25/64*x[3]**8)/(1-x[3]**2)**7.5)*(1+self.ADt)))*31536000
#        Eqs[2]= (x[2]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*(((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)/(1-x[3]**2)**6)*(x[0]*np.cos(x[4])+self.ADt*x[1])-((1+31/2*x[3]**2+255/8*x[3]**4+185/16*x[3]**6+25/64*x[3]**8)/(1-x[3]**2)**7.5)*(1+self.ADt)))*31536000
#        Eqs[3]= (x[3]*(27*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*((((1+3/2*x[3]**2+1/8*x[3]**4)/(1-x[3]**2)**5)*11/18*(x[0]*np.cos(x[4])+self.ADt*x[1])-((1+15/4*x[3]**2+15/8*x[3]**4+5/64*x[3]**6)/(1-x[3]**2)**6.5)*(1+self.ADt))))*31536000   
#        Eqs[4]= -3/2*(TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.mu()/self.Ms)**(1/2)*(x[0]*np.cos(2)+self.ADt*x[1])*np.sin(x[4])*((1+3*x[3]**2+3/8*x[3]**4)/(1-x[3]**2)**(5))*31536000
        
        return Eqs
    def i_Dt_model(self,t,x):
        Eqs= np.zeros((5))
        Eqs[0]= (-3*TenTai.G/(self.Cp()*x[2]**6)*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-1)*(((1+3*x[3]**2+3/8*x[3]**4)/(1-x[3]**2)**(9/2))*x[0]*np.cos(x[4])-((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)/(1-x[3]**2)**6))
                +3/2*x[0]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*(((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)/(1-x[3]**2)**6)*(x[0]*np.cos(x[4])+self.ADt*x[1]*np.cos(x[4]))-((1+31/2*x[3]**2+255/8*x[3]**4+185/16*x[3]**6+25/64*x[3]**8)/(1-x[3]**2)**7.5)*(1+self.ADt)))*31536000
        Eqs[1]= (-3*TenTai.G/(self.Cs()*x[2]**6)*self.k2s*self.Dts()*self.Mp**2*(self.Rs**5/self.Rp**6)*(((1+3*x[3]**2+3/8*x[3]**4)/(1-x[3]**2)**(9/2))*x[1]*np.cos(x[4])-((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)/(1-x[3]**2)**6))
                +3/2*x[1]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*(((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)/(1-x[3]**2)**6)*(x[0]*np.cos(x[4])+self.ADt*x[1]*np.cos(x[4]))-((1+31/2*x[3]**2+255/8*x[3]**4+185/16*x[3]**6+25/64*x[3]**8)/(1-x[3]**2)**7.5)*(1+self.ADt)))*31536000
        Eqs[2]= (x[2]*(6*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*(((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)/(1-x[3]**2)**6)*(x[0]*np.cos(x[4])+self.ADt*x[1]*np.cos(x[4]))-((1+31/2*x[3]**2+255/8*x[3]**4+185/16*x[3]**6+25/64*x[3]**8)/(1-x[3]**2)**7.5)*(1+self.ADt)))*31536000
        Eqs[3]= (x[3]*(27*TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-3)*((((1+3/2*x[3]**2+1/8*x[3]**4)/(1-x[3]**2)**5)*11/18*(x[0]*np.cos(x[4])+self.ADt*x[1]*np.cos(x[4]))-((1+15/4*x[3]**2+15/8*x[3]**4+5/64*x[3]**6)/(1-x[3]**2)**6.5)*(1+self.ADt))))*31536000   
        Eqs[4]= (-3*TenTai.G/(self.Cp()*x[2]**6)*self.k2p*self.Dtp*self.Ms**2*(self.Rp**-1)*(1-x[3]**2)**-6*1/x[0]*np.sin(x[4])*((1+15/2*x[3]**2+45/8*x[3]**4+5/16*x[3]**6)-1/2*(np.cos(x[4])-self.Cp()/self.mu()*x[2]**-2*(1-x[3]**2)**-0.5*x[0])*x[0]*(1-x[3]**2)**(3/2)*(1+3*x[3]**2+3/8*x[3]**4)))*31536000 
                
        #Eqs[4]= -3/2*(TenTai.G/(self.mu()*x[2]**8))*self.k2p*self.Dtp*self.Ms**2*(self.mu()/self.Ms)**(1/2)*(x[0]*np.cos(2)+self.ADt*x[1])*np.sin(x[4])*((1+3*x[3]**2+3/8*x[3]**4)/(1-x[3]**2)**(5))*31536000
        
        return Eqs








        