#
FC = FC
#  Check array bounds at run time:
FFLAGS =  -g -C 
FFLAGS = -fast -Wl,-Bstatic
FFLAGS = 
 
.f: 
	$(FC) $(FFLAGS) $*.f  

obj6 =	mdljmix4.o
mdljmix4  : $(obj6)
	$(FC) $(FFLAGS) $(LFLAGS) $(obj6) -o $@.x

 
clean :
	rm -f *.o
