library(pracma) # for evaluating cross product
library(OceanView) # for plotting curvature vectors

# growth function
rate<-function(r,K,x){
  return(r*x*(1-x/K))
}

# setting up environmental gradient; when b=0 we get a linear niche limit and when b!=0 we get a curved niche limit
environment <- function(b,N){
  x_start = 6; x_end = -6
  del_x = (x_end - x_start)/(N+1)
  y_start = -2.7; y_end = 2.7
  del_y = (y_end - y_start)/(N+1)

  env = matrix(0,N+2,N+2)
  for (i in 1:(N+2)){
    for (j in 1:(N+2)){
      env[j,i] = x_start + del_x*(i-1) + b*exp(-3*(y_start -1+ del_y*(j-1))^2) - b*exp(-3*(y_start+ 1+ del_y*(j-1))^2)}}
  return(env)
}


# Features of the landscape
dh = 1; # inter-lattice distance
length_x = 150; # Total x distance
n = floor(length_x/dh); # lattice rows/columns
x=(1:n)*dh;
y=(1:n)*dh;

# environment
b=3
alpha = .1; # change in growth rate per unit change in environment (aka sensitivity towards environment)
e_star = 5 # niche limit. If e<e_star species have a positive per-capita growth rate
e = environment(b,n) + e_star # environmental conditions on the geographical landscape
G_rate = alpha*(e-e_star) # per-capita growth rate on the geographical landscape

# model and simulation parameters
K=10; # carrying capacity
toll = 1e-5; # tolerance level to truncate the simulation
N_star = K/2; # population density at range limit (half-carrying capacity)
D = 25; # diffusion constant
dt = min(0.8*(dh*dh)/(4*D),0.01); # time step (denominator is 4 in 2D)


past_density = matrix(0,n+2,n+2)
curr_density = matrix(0,n+2,n+2)
curr_density[which(G_rate>0)]=K; # initializing regions with positive per-capita growth rate at carrying capacity and zero otherwise

minX = curr_density
plusX = curr_density
minY = curr_density
plusY = curr_density
I = 2:n+1; # central indices in y direction
J = 2:n+1; # central indices in x direction

# simulating Fisher equation
while (sum(abs(curr_density-past_density))>toll){
  past_density=curr_density;

  minX[I,J]=past_density[I,J-1];
  plusX[I,J]=past_density[I,J+1];
  minY[I,J]=past_density[I-1,J];
  plusY[I,J]=past_density[I+1,J];
  
  reaction = rate(G_rate,K,curr_density) #reaction term in Eq.1 of the paper
  diffusion = (D/(dh^2))*(plusX+minX+plusY+minY-4*past_density) #diffusion term in Eq.1 of the paper
  # total change in density
  curr_density = curr_density + dt*(reaction + diffusion);


  # Reflective boundary condition
  curr_density[,1] = curr_density[,2];
  curr_density[,n+2]= curr_density[,n+1];
  curr_density[1,]= curr_density[2,];
  curr_density[n+2,]= curr_density[n+1,];
}
I = (2*n/10):(8*n/10); # we leave out 20% of cells on the boundary of the lattice to remove any boundary effects due to reflective boundary conditions
J = I;

curr_density=t(curr_density[I,J])
G_rate=t(G_rate[I,J])
e = t(e[I,J])

Boundary=contourLines(1:length(I),1:length(I),curr_density, levels=c(N_star))
Boundary = list(x = Boundary[[1]]$x,y = Boundary[[1]]$y)# x and y coordinates of the boundary
BioClim=contourLines(1:length(I),1:length(I),G_rate, levels=c(0))
BioClim = list(x = BioClim[[1]]$x,y = BioClim[[1]]$y)# x and y coordinates of the niche limit, where species per-capita growth rate is zero

Boundary = as.data.frame(Boundary)

# Calculating curvature at the niche limit
dsx = diff(BioClim$x);
dsy = diff(BioClim$y);
ds = sqrt(dsx^2+dsy^2);
Tx = dsx/ds;
Ty = dsy/ds;

ds2 = 0.5*(ds[c(length(ds),1:(length(ds)-1))]+ds);
Hx = diff(Tx[c(length(ds),1:length(ds))])/ds2;
Hy = diff(Ty[c(length(ds),1:length(ds))])/ds2;

tangent = cbind(Tx,Ty,array(0,length(Tx)));
normal =  cbind(Hx,Hy,array(0,length(Hx)));
curvature = -cross(tangent,normal);
curvature = curvature[,3]/dh;
sign = array(0,length(curvature))
sign[which(curvature>0)] = curvature[which(curvature>0)]/max(curvature)
sign[which(curvature<0)] = -curvature[which(curvature<0)]/min(curvature)


# plotting figure 2 of the paper
par(mfrow=c(1,2),mai=c(0.35,0.35,0.35,0.02))

# when niche limit is a straight line
I = (2*n/10):(8*n/10); # central indices in y direction
J = I;

plot(0,0,xlim=c(5,length(I)*dh-5),ylim=c(5,length(I)*dh-5),cex=0,
     xlab="", ylab="",xaxt='n',yaxt='n',asp=1)
B = as.data.frame(rbind(c(1,1),c(1,length(I)),c(floor(length(I)/2),length(I)),c(floor(length(I)/2),1)))
colnames(B)<-c("x","y")
polygon(B$x*dh, B$y*dh,col="grey",border="grey")
points(c(floor(length(I)/2),floor(length(I)/2)), c(1,length(I)),type="l",col="black",lwd=1) # niche limit
mtext(side = 1, "Space(x)",cex = 1,line = 0.5)
mtext(side = 3, "Occupied\nregion",cex = 0.8, at=length(I)/4 ,line = -3)
mtext(side = 3, "Unoccupied\nregion",cex = 0.8, at=(3/4)*length(I) ,line = -3)
mtext(side = 2, "Space(y)",cex = 1,line = 0.5)
mtext(side = 3, "Linear Geometry",cex = 1,line = 0.5)
mtext(side = 3, "A",cex = 1.2, at =-5, line = 0.5)

# when niche limit is a curved contour
plot(0,0,xlim=c(5,length(I)*dh-5),ylim=c(5,length(I)*dh-5),cex=0,
     xlab="", ylab="",xaxt='n',yaxt='n',asp=1)
B = rbind(c(1,Boundary$y[1]),Boundary,c(1,Boundary$y[length(Boundary$y)]))
polygon(B$x[1:length(ds)]*dh, B$y[1:length(ds)]*dh,col="grey",border="grey")
arrow = as.data.frame(list(Hx = 60*Hx,Hy = 60*Hy,x = BioClim$x[1:length(ds)]*dh, y =BioClim$y[1:length(ds)]*dh,sign = sign))
arrow$mod = sqrt(arrow$Hx^2+arrow$Hy^2)
arrow = arrow[-c(1,which(arrow$x<67 & arrow$x>28)),]
arrow = arrow[-c(1,which(arrow$mod<1)),]
arrow = rbind(arrow,c(-10, 0, 41.5, 18.5, 0,20))
vectorplot(arrow$Hx,arrow$Hy,arrow$x, arrow$y,by=1,arr=T,arr.width=.05,arr.length = 0.07,add=T,lwd=0.75,col="grey30") # plotting curvature vector
points(BioClim$x*dh, BioClim$y*dh,type="l",col="black",lwd=1) # niche limit
mtext(side = 1, "Space(x)",cex = 1,line = 0.5)
mtext(side = 3, "Curved Geometry",cex = 1,line = 0.5)
mtext(side = 3, "B",cex = 1.2, at =-5, line = 0.5)
legend(0.3*floor(length(I)),0.7*floor(length(I)/2), legend=c("Niche limit", "Curvature vector"), col=c("black","transparent"), bg="transparent",lty=c(1,1), cex=0.75, text.font=0.9,box.lty=0)


# Calculating curvature at the species boundary
dsx = diff(Boundary$x);
dsy = diff(Boundary$y);
ds = sqrt(dsx^2+dsy^2);
Tx = dsx/ds;
Ty = dsy/ds;

ds2 = 0.5*(ds[c(length(ds),1:(length(ds)-1))]+ds);
Hx = diff(Tx[c(length(ds),1:length(ds))])/ds2;
Hy = diff(Ty[c(length(ds),1:length(ds))])/ds2;

tangent = cbind(Tx,Ty,array(0,length(Tx)));
normal =  cbind(Hx,Hy,array(0,length(Hx)));
curvature = -cross(tangent,normal);
curvature = curvature[,3]/dh;
sign = array(0,length(curvature))
sign[which(curvature>0)] = curvature[which(curvature>0)]/max(curvature)
sign[which(curvature<0)] = -curvature[which(curvature<0)]/min(curvature)

# plotting figure S1 of the paper
par(mfrow=c(1,1),mai=c(0.7,.7,0.35,0.05))
index_boundary = round(Boundary$x)+(round(Boundary$y)-1)*91
index_boundary = index_boundary[1:(length(index_boundary)-1)]
dat = list(env =e[index_boundary], curv = abs(curvature))
plot(dat$env, dat$curv, xlab="", ylab="", pch = 19,cex= 0.5)
mtext(side = 1, "Environment",cex = 1,line = 2)
mtext(side = 2, "Absolute Curvature",cex = 1,line = 2)
