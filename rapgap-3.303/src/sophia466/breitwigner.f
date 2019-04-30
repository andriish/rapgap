

       DOUBLE PRECISION function breitwigner(sigma_0,Gamma,
     &                     DMM,eps_prime)

       IMPLICIT DOUBLE PRECISION (A-M,O-Z)
       IMPLICIT INTEGER (N)

       SAVE

c***************************************************************************
c calculates Breit-Wigner cross section of a resonance with width Gamma [GeV],
c mass DMM [GeV], max. cross section sigma_0 [mubarn] and total mass of the 
c interaction s [GeV] 
c***************************************************************************
       pm = 0.93827D0
       s = pm*pm+2.D0*pm*eps_prime
       gam2s = Gamma*Gamma*s
       breitwigner = sigma_0
     &              *(s/eps_prime**2)*gam2s/((s-DMM*DMM)**2+gam2s)

       RETURN
       
       END
