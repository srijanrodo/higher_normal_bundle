import random
import time

F = QQ
R.<x,y> = F[]
  
def check_for_this_curve(curve, k=2):
    r = len(curve)-1
    d = curve[0].degree()
    derivatives = matrix([[f.derivative(x,i,y,k-i) for i in range(k+1)] for f in curve])

    #Check for degenerate behavior
    minors = derivatives.minors(k+1)
    I = R.ideal(minors)
    J = R.ideal(x,y)
    if not I.radical() == J:
        print('The curve has the wrong osculating behaviour')
        return(-1)

    p = ((k+1)*(d-k))//(r-k) #Expected small degree of the relation 
    a = (r-k) - (((k+1)*(d-k))%(r-k)) #Expected number of small degree relations

    def compute_defect(derivatives,s):
        Q = []
        for l in range(k+1):
            for m in range(d-k+s+1):
                row = []
                for j in range(r+1):
                    for i in range(s+1):
                        row.append(F(derivatives[j][l].coefficient({x:m-i,y:d-k-m+i})))
                Q.append(list(row))
        Q = matrix(Q)
        K = Q.right_kernel().basis()
        
        relations=[]
        for v in K:
            l=[]
            index=0
            for j in range(r+1):
                f = 0
                for i in range(s+1):
                    f = f + v[index]*x^i*y^(s-i)
                    index=index+1
                l.append(f)
            relations.append(list(l))
        #print(relations)
        return (Q.ncols() - Q.rank(), list(relations))

    if not compute_defect(derivatives,p-1)[0] == 0:
        print('There is a relation in degree '+str(p-1)+', which should not happen')
        return(-2)
    res=compute_defect(derivatives,p)
    if not res[0] == a:
        print('There are too many relations in degree',p)
        return(-3)
    print('BALANCED!')
    return (1,res[1])

def generate_rand_equations(d,seed=0,rang=200,size=-1):
    if size == -1:
        size=d+1
    random.seed(seed)
    f = 0
    indices=list(range(d+1))
    random.shuffle(indices)
    for i in indices:
        f = f + random.randint(0,rang)*(x^i*y^(d-i))
        if len(f.terms()) >= size:
            break
    return f
#curve = [x^7*y + y^8, x^4*y^4 + x^2*y^6, x*y^7, x^8, x^6*y^2]
#check_for_this_curve(curve)
r=4
for deg in range(2001,2002):
    with open('balanced_ab1','a') as f:
        print("For degree ",deg,file=f)
    i = (deg // 2) - 1
    curve=[x^i*y^(deg-i), x^(deg-2)*y^2 + x^2*y^(deg-2), x^deg,x^(deg-1)*y+x*y^(deg-1), y^(deg)]
    res=check_for_this_curve(curve)
    #print('hi', res)
    if res[0] == 1:
        with open('balanced_ab1','a') as f:
            #print(curve[0],file=f)
            for i in res[1][0]:
                print(i,file=f)
            print('\n',file=f)
    else:
       assert False
'''
for i in range(100):
    curve=[ generate_rand_equations(deg, seed=time.time(),rang=100) for j in range(r+1)]#+ [x^(deg-2)*y^2 + x^2*y^(deg-2),x^deg,x^(deg-1)*y+x*y^(deg-1),y^(deg)] #[generate_rand_equations(deg, seed=time.time(),rang=1,size=1) for j in range(r-1)]
    res=check_for_this_curve(curve)
    if res == 1:
        with open('balanced','a') as f:
            print(curve,'\n',file=f)
        #deg=deg+1
        #continue
    elif res == -1:
        with open('degenerate','a') as f:
            print(curve,'\n',file=f)
    elif res == -2:
        with open('toolow','a') as f:
            print(curve,'\n',file=f)
    elif res == -3:
        with open('toohigh','a') as f:
            print(curve,'\n',file=f)
'''