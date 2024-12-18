#procedure Kinematics

id p1.p1 = 0;
id p2.p2 = 0;
id p3.p3 = 0;
id q.q = 0;

#do i=1,3
  id p`i'.p`i' = 0;
  id p`i'.q = s`i'q/2;
  #do j=`i'+1,3
    id p`i'.p`j' = s`i'`j'/2;
  #enddo
#enddo

id q.E = 0;

.sort

#endprocedure