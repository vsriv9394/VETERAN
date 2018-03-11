X   =   [i for i=0:0.01:10,j=0:1:1,k=0:1]
Y   =   [j for i=0:0.01:10,j=0:1:1,k=0:1]
Z   =   [k for i=0:0.01:10,j=0:1:1,k=0:1]

num =   size(X,1)*size(X,2)*size(X,3)

Xn  =   reshape(X,num)
Yn  =   reshape(Y,num)
Zn  =   reshape(Z,num)

Xnum    =   size(X,1)
Ynum    =   size(X,2)
Znum    =   size(X,3)

Xn  =   vcat(Xnum,Xn)
Yn  =   vcat(Ynum,Yn)
Zn  =   vcat(Znum,Zn)

writedlm("new.mesh",hcat(Xn,Yn,Zn))

BCArr = zeros(2*(Ynum[1]-1)*(Znum[1]-1)+2*(Znum[1]-1)*(Xnum[1]-1)+2*(Xnum[1]-1)*(Ynum[1]-1),4)
m     = 0
for j=1:Ynum[1]-1
    for k=1:Znum[1]-1
        m       =   m+1
        BCArr[m,1]  =   1
        BCArr[m,2]  =   j
        BCArr[m,3]  =   k
        BCArr[m,4]  =   -1
    end
end
for j=1:Ynum[1]-1
    for k=1:Znum[1]-1
        m       =   m+1
        BCArr[m,1]  =   Xnum[1]
        BCArr[m,2]  =   j
        BCArr[m,3]  =   k
        BCArr[m,4]  =   -1
    end
end

for i=1:Xnum[1]-1
    for k=1:Znum[1]-1
        m       =   m+1
        BCArr[m,1]  =   i
        BCArr[m,2]  =   1
        BCArr[m,3]  =   k
        BCArr[m,4]  =   -1
    end
end
for i=1:Xnum[1]-1
    for k=1:Znum[1]-1
        m       =   m+1
        BCArr[m,1]  =   i
        BCArr[m,2]  =   Ynum[1]
        BCArr[m,3]  =   k
        BCArr[m,4]  =   -1
    end
end

for i=1:Xnum[1]-1
    for j=1:Ynum[1]-1
        m       =   m+1
        BCArr[m,1]  =   i
        BCArr[m,2]  =   j
        BCArr[m,3]  =   1
        BCArr[m,4]  =   -1
    end
end
for i=1:Xnum[1]-1
    for j=1:Ynum[1]-1
        m       =   m+1
        BCArr[m,1]  =   i
        BCArr[m,2]  =   j
        BCArr[m,3]  =   Znum[1]
        BCArr[m,4]  =   -1
    end
end

writedlm("new.bc",BCArr)

Cells = zeros(size(X,1)+1,size(X,2)+1,size(X,3)+1,6)
Cells[:,:,:,2] = zeros(size(X,1)+1,size(X,2)+1,size(X,3)+1,1)
Cells[:,:,:,3] = zeros(size(X,1)+1,size(X,2)+1,size(X,3)+1,1)
Cells[:,:,:,4] = zeros(size(X,1)+1,size(X,2)+1,size(X,3)+1,1)

ind                             = convert(Int,floor(size(X,1)/2.0))

Cells[1:ind,:,:,5]              = ones(ind,size(X,2)+1,size(X,3)+1,1)
Cells[ind+1:size(X,1)+1,:,:,5]  = 0.1*ones(size(X,1)+1-ind,size(X,2)+1,size(X,3)+1,1)

Cells[1:ind,:,:,1]              = ones(ind,size(X,2)+1,size(X,3)+1,1)
Cells[ind+1:size(X,1)+1,:,:,1]  = 0.125*ones(size(X,1)+1-ind,size(X,2)+1,size(X,3)+1,1)

Cells[:,:,:,6]                  = ones(size(X,1)+1,size(X,2)+1,size(X,3)+1,1)
write("initcond.vtrn",Cells)
