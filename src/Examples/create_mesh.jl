X   =   [i for i=0:0.01:10,j=0:0.01:1,k=0:0.01:0.01]
Y   =   [j for i=0:0.01:10,j=0:0.01:1,k=0:0.01:0.01]
Z   =   [k for i=0:0.01:10,j=0:0.01:1,k=0:0.01:0.01]

num =   1001*101*2

Xn  =   reshape(X,num)
Yn  =   reshape(Y,num)
Zn  =   reshape(Z,num)

Xnum    =   [1001,]
Ynum    =   [101,]
Znum    =   [2,]

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
        BCArr[m,4]  =   -2
    end
end

for i=1:Xnum[1]-1
    for k=1:Znum[1]-1
        m       =   m+1
        BCArr[m,1]  =   i
        BCArr[m,2]  =   1
        BCArr[m,3]  =   k
        BCArr[m,4]  =   -3
    end
end
for i=1:Xnum[1]-1
    for k=1:Znum[1]-1
        m       =   m+1
        BCArr[m,1]  =   i
        BCArr[m,2]  =   Ynum[1]
        BCArr[m,3]  =   k
        BCArr[m,4]  =   -3
    end
end

for i=1:Xnum[1]-1
    for j=1:Ynum[1]-1
        m       =   m+1
        BCArr[m,1]  =   i
        BCArr[m,2]  =   j
        BCArr[m,3]  =   1
        BCArr[m,4]  =   -3
    end
end
for i=1:Xnum[1]-1
    for j=1:Ynum[1]-1
        m       =   m+1
        BCArr[m,1]  =   i
        BCArr[m,2]  =   j
        BCArr[m,3]  =   Znum[1]
        BCArr[m,4]  =   -3
    end
end

writedlm("new.bc",BCArr)
