using LinearAlgebra,Printf,MKL
const m0=1.0::Float64
const t=1.0::Float64
const μ=0.0::Float64
const Δ=0.4::Float64
const λ=1.0::Float64
# const J=1.0::Float64
const Bz=0.0::Float64
const g=1.382301::Float64
const ORB=8::Int64
const N1=50::Int64
const N2=ORB*N1::Int64
const N=(ORB*N1*N1)::Int64
const Nf=N1*N1::Int64
const Name="OBC"::String
BLAS.set_num_threads(10)
function BuiltH(H,J)
    σx=[0 1;1 0]
    σy=[0 -im;im 0]
    σz=[1 0;0 -1]
    σ0=[1 0;0 1]
    # PH (x) Orbital (x) Spin
    Γ1=kron(σz,σx,σz) #SOC x
    Γ2=kron(σz,σy,σ0) #SOC y
    Γ3=kron(σz,σz,σ0) #Hopping
    Γ4=kron(σx,σ0,σ0) #SC
    Γ6=kron(σz,σ0,σ0) #Chemical Pot
    Γx=kron(σ0,σ0,σx) 
    Γy=kron(σ0,σ0,σy) 
    Γz=kron(σ0,σ0,σz) 
    σx=0
    σy=0
    σz=0

    σ0=0

    

    # data = readdlm("BlochSkyrmion_SpinConfig.txt", '\t')
    # data = readdlm("NeelSkyrmion_SpinConfig.txt", '\t')

    # # Extract Spin_x, Spin_y, Spin_z
    # Spin_x = data[:, 4]
    # Spin_y = data[:, 5]
    # Spin_z = data[:, 6]

    # # Reshape arrays
    # Spin_x = reshape(Spin_x, N1, N1)
    # Spin_y = reshape(Spin_y, N1, N1)
    # Spin_z = reshape(Spin_z, N1, N1)


    
    # for l in 0:N1-1,j in 1:ORB:N2
    #     y=l+1
    #     x=(j÷ORB)+1
    #     for col in 1:ORB,row in 1:ORB
    #         H[((l*N2)+j+row-1),((l*N2)+j+col-1)]=μ*Γ1[row,col]+Δ*Γ4[row,col]+J*sin(g*x)*Γx[row,col]+J*sin(g*y)*Γy[row,col]+(J*cos(g*x)+J*cos(g*y)+B)*Γz[row,col]
    #         # H[((l*N2)+j+row-1),((l*N2)+j+col-1)]=μ*Γ1[row,col]+Δ*Γ4[row,col]+J*Γx[row,col]
    #         # H[((l*N2)+j+row-1),((l*N2)+j+col-1)]=μ*Γ1[row,col]+Δ*Γ4[row,col]+J*cos((x+y)*g)*Γx[row,col]+J*sin((x+y)*g)*Γy[row,col]
    #         # H[((l*N2)+j+row-1),((l*N2)+j+col-1)]=μ*Γ1[row,col]+Δ*Γ4[row,col]+J*Spin_x[x,y]*Γx[row,col]+J*Spin_y[x,y]*Γy[row,col]+J*Spin_z[x,y]*Γz[row,col]

    #         if j+(ORB+1)<=ORB*N1
    #             H[((l*N2)+j+row-1),((l*N2)+j+col-1+ORB)]=t*Γ1[row,col]+Δd*Γ4[row,col]
    #         end
    #         if l<N1-1
    #             H[((l*N2)+j+row-1),(((l+1)*N2)+j+col-1)]=t*Γ1[row,col]-Δd*Γ4[row,col]
    #         end
    #     end
    # end
    for l in 0:N1-1,j in 1:ORB:N2
        My=l+1
        Mx=(j÷ORB)+1
        for col in 1:ORB,row in 1:ORB
            H[((l*N2)+j+row-1),((l*N2)+j+col-1)]=m0*Γ3[row,col]+Δ*Γ4[row,col]-μ*Γ6[row,col]+J*sin(g*My)*Γy[row,col]+(J*cos(g*Mx)+J*cos(g*My)+Bz)*Γz[row,col]
            if j+(ORB+1)<=ORB*N1
                H[((l*N2)+j+row-1),((l*N2)+j+col-1+ORB)]=-(im*λ)*Γ1[row,col]-t*Γ3[row,col]
            end
            if l<N1-1
                H[((l*N2)+j+row-1),(((l+1)*N2)+j+col-1)]=-(im*λ)*Γ2[row,col]-t*Γ3[row,col]
            end
        end
    end
    # for i in 0:Nf-1
    #     x=i%N1
    #     y=i÷N1
    #     xn=(x+1)%N1
    #     yn=(y+1)%N1
    #     xp=y*N2+x*ORB
    #     xnp=(y*N2)+xn*ORB
    #     yp=y*N2+x*ORB
    #     ynp=yn*N2+x*ORB
    #     Mx=x+1
    #     My=y+1
    #     #===============OBC===============#
    #     for col in 1:ORB,row in 1:ORB 
    #         #(x,y)→(x,y) :onsite
    #         H[xp+row,xp+col]=m0*Γ3[row,col]+Δ*Γ4[row,col]-μ*Γ6[row,col]+J*sin(g*My)*Γy[row,col]+(J*cos(g*Mx)+J*cos(g*My)+Bz)*Γz[row,col]
    #         #(x,y)→(x+1,y) :hopping along x
    #         if xn>x
    #             H[xp+row,xnp+col]=-(im*λ)*Γ1[row,col]-t*Γ3[row,col]
    #             # H[xnp+col,xp+row]=conj(-(im*λ)*Γ1[row,col]-t*Γ3[row,col])
    #         end
    #         #(x,y)→(x,y+1) :hopping along y
    #         if yn>y
    #             H[xp+row,ynp+col]=-(im*λ)*Γ2[row,col]-t*Γ3[row,col]
    #             # H[ynp+col,xp+row]=conj(-(im*λ)*Γ2[row,col]-t*Γ3[row,col])
    #         end
    #     end
       
        
    # end
    Γ1=0
    Γ2=0
    Γ3=0
    Γ4=0
end
function doscalculate(lState,H,DOS)
    # δ=0.01
    l=lState
    for i in 0:N1-1,j in 0:N1-1
        tr=0.0
        for p in (ORB*i+j*N2+1):(ORB*i+j*N2+ORB)
            tr+=H.vectors[p,l]*conj(H.vectors[p,l])
        end
        # DOS[j+1,i+1]+=abs(tr*δ)/(((E-nu[l])*(E-nu[l]))+(δ*δ));
        DOS[j+1,i+1]+=abs(tr)
    end
end
function main()
    H=zeros(ComplexF64,N,N)
    fp1=open("eigenvalueJ$Name.txt","w")
    parameter=1.440000  
    BuiltH(H,parameter)
    H=eigen(Hermitian(H))
    for i in 1:N 
        @printf(fp1,"%d %f \n",i,H.values[i])
    end
    fp1=open("dos$Name.txt","w")
    DOS=zeros(Float64,N1,N1)
    del=1e-3
    E=0.0
    ls=1
    le=N
    # for l in ls:le,k in 0:N1-1,i in 0:N1-1
    #     tr=0.0;
    #     for p in (ORB*i+k*N2+1):(ORB*i+k*N2+ORB)
    #         tr+=H.vectors[p,l]*conj(H.vectors[p,l]);
    #     end
    #     # DOS[k+1,i+1]+=abs(tr)
    #     DOS[k+1,i+1]+=(abs(tr)*del)/((E-real(H.values[l]))^2+del^2);
    # end
    for l in ls:le,i in 0:Nf-1
        x=i%N1
        y=i÷N1
        xp=y*N2+x*ORB
        tr=0.0
        for p in xp+1:xp+ORB
            tr+=H.vectors[p,l]*conj(H.vectors[p,l])
        end
        DOS[x+1,y+1]+=(abs(tr)*del)/((E-real(H.values[l]))^2+del^2)
    end
    MaxDOS=maximum(DOS)::Float64
    for i in 1:N1,j in 1:N1
        @printf(fp1,"%d %d %f \n",j,i,DOS[i,j]/MaxDOS)
    end

    DOS=Nothing
    close(fp1)
 
    
    for lState in (N÷2)-3:(N÷2)+4
        fp1=open("dosState$lState.txt","w")
        DOS=zeros(Float64,N1,N1)
        doscalculate(lState,H,DOS)
        MaxDOS=maximum(DOS)::Float64
        for j in 1:N1,i in 1:N1
            @printf(fp1,"%d %d %f \n",i,j,DOS[i,j]/MaxDOS)
        end
        DOS=0
        close(fp1)
    end
    H=Nothing
    

end
main()
function mainLoop()
    
    fp1=open("eigenvalueLoop$Name.txt","w")
    for l in 0:50 
        parameter=0.0+2.0*l/50
        H=zeros(ComplexF64,N,N)
        BuiltH(H,parameter)
        H=eigvals(Hermitian(H))
        for i in 1:N 
            @printf(fp1,"%f %f \n",parameter,H[i])
        end
        H=0
    end
    
    
    close(fp1)
end
# mainLoop()