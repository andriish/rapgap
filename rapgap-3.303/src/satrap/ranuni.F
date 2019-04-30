      FUNCTION RANUNI(A,B)
C ... Uniform random random number in range [A,B]
      DOUBLE PRECISION RANUNI,A,B,RN
	Double Precision draprn

ccc      Call RANMAR(RVec,1)
c      RN=RVec(1)
c changed to standard draprn

      RN=draprn()
      RANUNI=A+RN*(B-A)
      END
