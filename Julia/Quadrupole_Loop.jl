using LinearAlgebra,Printf
const m0=1.0::Float64
const t=1.0::Float64
const λ=1.0::Float64
const μ=0.0::Float64
const Δ=0.4::Float64
const Bz=0.5::Float64
# const J=0.8::Float64
# const g=π/4.0::Float64
const ORB=8::Int64
const N1=16::Int64
const N2=ORB*N1::Int64
const N=(ORB*N1*N1)::Int64
const M=(N÷2)::Int64
import Base.Threads.@threads
BLAS.set_num_threads(1)
function BuiltHam(H,J,g)
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
        #PBC along x
        # if j==1
        #     H[((l*N2)+row),((l*N2)+(N1-1)*ORB+col)]=(im*λ)*Γ1[row,col]-(im*t*g/2.0)*Γ7[row,col]-t*Γ3[row,col]
        # end
        #PBC along y
        # if l==0
        #     H[(j+row-1),(((N1-1)*N2)+j+col-1)]=(im*λ)*Γ2[row,col]-(im*t*g/2.0)*Γ7[row,col]-t*Γ3[row,col]
        # end
    end
    Γ1=0
    Γ2=0
    Γ3=0
    Γ4=0
end
function Quadrupole(parameter1,parameter2)::Float64
    H=zeros(ComplexF64,N,N)
    BuiltHam(H,parameter1,parameter2)
    H=eigen(Hermitian(H))
    U=zeros(ComplexF64,N,M)
    for i in 1:M 
        U[:,i]=H.vectors[:,i]
    end
    H=0
    q=zeros(ComplexF64,N,N)
    for i in 0:N1-1,j in 0:N1-1
        for p in ORB*i+j*N2+1:ORB*i+j*N2+ORB
            q[p,p]=(i+1.0)*(j+1.0)/(N1*N1)
        end
    end 
    Q=exp(im*2.0*π*q)
    V=U'*Q*U
    U=0
    Q=0
    c1=det(V)
    c2=exp(-im*π*tr(q))
    qxy=(1.0/(2.0*π))*imag(log(c1*c2))
    Qxy=real(qxy)
    return Qxy
    V=0
end
function main()
    lN=50
    LL=(lN+1)*(lN+1)-1
    @threads for ll in 0:LL
        lx=ll%(lN+1)
        ly=ll÷(lN+1)
        parameter1=6.0*lx/lN  #J
        parameter2=π*ly/lN  #g
        Qxy=Quadrupole(parameter1,parameter2)
        fp=open("inv$ll.txt","w")
        @printf(fp,"%f %f %f %f \n",parameter1,parameter2,Qxy,abs(Qxy))    
        close(fp)
    end

end
main()
