# Paper Examples
from rnr import qnr
from numpy import array, imag, real, round
from matplotlib import pyplot as plt
from networkx import DiGraph, draw_shell

###############################################
###             save fig                    ###
###############################################
def save_fig(l,name,ratio,nrows,ncols):
    f,e = qnr(l)
    g = DiGraph(l)
    dx = 0.25; dy = 0.25
    
    fig, axes = plt.subplots(nrows=nrows,ncols=ncols)
    fig.tight_layout()
    
    lowx = int(round(min(real(e))-dx))
    highx = int(round(max(real(e))+dx))
    lowy = int(round(min(imag(e))-dy))
    highy = int(round(max(imag(e))+dy))
    
    axes[1].spines['bottom'].set_bounds((lowx,highx))
    axes[1].spines['left'].set_bounds((lowy,highy))
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    axes[1].plot(real(f),imag(f),color='#000000',fillstyle='right')
    axes[1].fill(real(f),imag(f),color='#C0C0C0')
    axes[1].plot(real(e),imag(e),color='#606060',linestyle='none',marker='*',markersize=10.0)
    axes[1].set_aspect(ratio)
    axes[1].set_xticks(range(lowx,highx+1))
    axes[1].set_yticks(range(lowy,highy+1))
    
    axes[0].axis('off')
    draw_shell(g,node_color='#C0C0C0',arrowsize=15,arrowstyle='-|>',ax=axes[0])
    axes[0].axis('equal')
    
    fig.savefig(name,dpi=400)
###############################################
###             main                        ###
###############################################
def main():
    try:
        # case 1 (cycle)
        l = array([[1,-1,0,0,0,0],[0,1,-1,0,0,0],[0,0,1,-1,0,0],[0,0,0,1,-1,0],[0,0,0,0,1,-1],[-1,0,0,0,0,1]],dtype=float)
        save_fig(l,"eps/case1.eps",0.6,2,1)
        # case 2 (directed join)
        l = array([[4,0,-1,-1,-1,-1],[-1,4,0,-1,-1,-1],[0,-1,4,-1,-1,-1],[0,0,0,1,-1,0],[0,0,0,-1,1,0],[0,0,0,0,0,0]],dtype=float)
        save_fig(l,"eps/case2.eps",1.5,2,1)
        # case 3 (partial directed join)
        l = array([[3,-1,0,0,-1,-1],[-1,3,0,0,-1,-1],[0,0,1,0,-1,0],[0,0,0,1,0,-1],[0,0,0,-1,1,0],[0,0,-1,0,0,1]],dtype=float)
        save_fig(l,"eps/case3.eps",1,2,1)
        # circulant normal
        l = array([[2,-1,-1,0,0,0],[0,2,-1,-1,0,0],[0,0,2,-1,-1,0],[0,0,0,2,-1,-1],[-1,0,0,0,2,-1],[-1,-1,0,0,0,2]],dtype=float)
        save_fig(l,"eps/circ-normal.eps",0.3,1,2)
        # balanced, non-normal
        l = array([[1,0,-1,0,0,0],[0,1,-1,0,0,0],[0,0,2,-1,-1,0],[0,0,0,1,0,-1],[0,0,0,0,1,-1],[-1,-1,0,0,0,2]],dtype=float)
        save_fig(l,"eps/bal-non-normal.eps",0.7,1,2)
        # 4-cycle
        l = array([[1,-1,0,0],[0,1,-1,0],[0,0,1,-1],[-1,0,0,1]],dtype=float)
        save_fig(l,"eps/4-cycle.eps",0.6,2,1)
        # twin splitting
        l = array([[1,0,-1,0,0,0],[0,1,-1,0,0,0],[0,0,2,-1,-1,0],[0,0,0,1,0,-1],[0,0,0,0,1,-1],[-1,-1,0,0,0,2]],dtype=float)
        save_fig(l,"eps/twin-splitting.eps",0.7,2,1)
        # restore normality
        l = array([[2,0,-1,0,0,0,-1],[0,2,-1,0,0,0,-1],[0,0,2,-1,-1,0,0],[0,0,0,2,0,-1,-1],[0,0,0,0,2,-1,-1],[-1,-1,0,0,0,2,0],[-1,-1,0,-1,-1,0,4]],dtype=float)
        save_fig(l,"eps/nrml-rest.eps",1,2,1)
        # non-djoin restricted normal
        l = array([[2,-1,0,0,-1,0,0,0],[-1,2,0,0,0,-1,0,0],[0,-1,3,-1,-1,0,0,0],[-1,0,-1,3,0,-1,0,0],[-1,-1,0,-1,4,0,0,-1],[-1,-1,-1,0,0,4,-1,0],[-1,-1,-1,0,-1,-1,5,0],[-1,-1,0,-1,-1,-1,0,5]],dtype=float)
        save_fig(l,"eps/rnrml-non-djoin.eps",1,1,2)
    except Exception as e:
        print(e)
if __name__ == '__main__':
    main()