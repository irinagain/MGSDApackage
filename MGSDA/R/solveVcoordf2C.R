.solveVcoordf2C <-
function(D,W,V,lambda,p,r,eps,maxiter)
{
  niter=0;
  .C("solveVcoordf2", as.double(as.vector(D)),as.double(as.vector(W)),as.double(as.vector(V)),as.double(lambda),as.integer(p),as.integer(r),as.double(eps),as.integer(maxiter),as.integer(niter))
}
