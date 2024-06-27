#procedure Dirac

id g_(0, mu1?, mu2?, mu3?, mu4?, mu5?) = g(mu1) * g(mu2) * g(mu3) * g(mu4) * g(mu5);
id g_(0, mu1?, mu2?, mu3?, mu4?) = g(mu1) * g(mu2) * g(mu3) * g(mu4);
id g_(0, mu1?, mu2?, mu3?) = g(mu1) * g(mu2) * g(mu3);
id g_(0, mu1?, mu2?) = g(mu1) * g(mu2);
id g_(0, mu1?) = g(mu1);
.sort

repeat;
    id disorder g(mu1?) * g(mu2?) = 2 * d_(mu1, mu2) - g(mu2) * g(mu1);
    renumber 1;
endrepeat;
* permute p3 to the left
#do i = 1,5;
    id g(mu1?) * g(p3) = 2 * p3(mu1) - g(p3) * g(mu1);
#enddo;

* permute p4 to the right
#do i = 1,5;
    id g(p4)*g(mu1?) = 2 * p4(mu1) - g(mu1)*g(p4);
#enddo;

id g(p1?) * g(p1?) = p1.p1;
id g(mu1?) * g(mu1?) = d;


#endprocedure