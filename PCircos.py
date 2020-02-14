import matplotlib.pyplot as plt
import math
import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import getopt
import sys

opts,args=getopt.getopt(sys.argv[1:],"I:i:h:")

fig,ax = plt.subplots(1)
listchrname=[]
listchrlength=[]
listhudu=[]
listexplode=[]
outer_colors=[]
def openchrfile(filename):
    filechr=open(filename,'r',encoding="utf-8")
    contentchr=filechr.read()
    contschr=contentchr.split('\n')
    a=0
    for i in contschr:
        if contschr[a].find('chr')>=0:
            chrname=contschr[a].split('\t')[2]
            listchrname.append(chrname)
            listchrlength.append(contschr[a].split('\t')[5])
            listexplode.append(0.05)
        a=a+1
    filechr.close()

radi=1
widt=0.15
vals = np.array([[60., 32.], [37., 40.], [29., 10.]])
list1=[]
listgene=[]
listgenedensity=[]
def calc_angle(x1,y1,x2,y2): 
    angle=0
    dy= y2-y1
    dx= x2-x1
    if dx==0 and dy>0:
        angle = 90
    if dx==0 and dy<0:
        angle = 270
    if dy==0 and dx>0:
        angle = 0
    if dy==0 and dx<0:
        angle = 180
    if dx>0 and dy>0:
        angle = math.atan(dx/dy)*180/math.pi
    elif dx<0 and dy>0:
        angle = -(90+math.atan(dx/dy)*180/math.pi)
    elif dx<0 and dy<0:
        angle = 360 + math.atan(dx/dy)*180/math.pi
    elif dx>0 and dy<0:
        angle = -(90+math.atan(dx/dy)*180/math.pi)
    return angle

def calposbyangle(angle,sign):
    if sign==1 and t==1:
            if angle<=90 and angle>=0:
                y=(radi+linkradi-widt/2)*math.sin(math.radians(angle))
                x=((radi+linkradi-widt/2)*(radi+linkradi-widt/2)-y*y)**0.5
            if 90<angle<=180:
                y=(radi+linkradi-widt/2)*math.sin(math.radians(180-angle))
                x=-((radi+linkradi-widt/2)*(radi+linkradi-widt/2)-y*y)**0.5
            if 180<angle<270:
                y=-(radi+linkradi-widt/2)*math.sin(math.radians(angle-180))
                x=-((radi+linkradi-widt/2)*(radi+linkradi-widt/2)-y*y)**0.5
            if angle<=360 and angle>=270:
                y=-(radi+linkradi-widt/2)*math.sin(math.radians(360-angle))
                x=((radi+linkradi-widt/2)*(radi+linkradi-widt/2)-y*y)**0.5
            list1.append(x)
            list1.append(y)
    else:
        if angle<=90 and angle>=0:
            y=(radi-widt/2)*math.sin(math.radians(angle))
            x=((radi-widt/2)*(radi-widt/2)-y*y)**0.5
        if 90<angle<=180:
            y=(radi-widt/2)*math.sin(math.radians(180-angle))
            x=-((radi-widt/2)*(radi-widt/2)-y*y)**0.5
        if 180<angle<270:
            y=-(radi-widt/2)*math.sin(math.radians(angle-180))
            x=-((radi-widt/2)*(radi-widt/2)-y*y)**0.5
        if angle<=360 and angle>=270:
            y=-(radi-widt/2)*math.sin(math.radians(360-angle))
            x=((radi-widt/2)*(radi-widt/2)-y*y)**0.5
        if sign==1:
            list1.append(x)
            list1.append(y)
        if sign==0:
            ax.annotate('c', xy=(x, y),rotation=90-angle)
        if sign==2:
            listgene.append(x)
            listgene.append(y)
        if sign==3:
            listgenedensity.append(x)
            listgenedensity.append(y)


def paintchrbyangle(angle,string):
    radichr=1
    if angle<=90 and angle>=0:
        y=(radichr-0.15)*math.sin(math.radians(angle))
        x=((radichr-0.15)*(radichr-0.15)-y*y)**0.5
    if 90<angle<=180:
        y=(radichr-0.15)*math.sin(math.radians(180-angle))
        x=-((radichr-0.15)*(radichr-0.15)-y*y)**0.5
    if 180<angle<270:
        y=-(radichr-0.15)*math.sin(math.radians(angle-180))
        x=-((radichr-0.15)*(radichr-0.15)-y*y)**0.5
    if angle<=360 and angle>=270:
        y=-(radichr-0.15)*math.sin(math.radians(360-angle))
        x=((radichr-0.15)*(radichr-0.15)-y*y)**0.5
    rets=str(x)+'\t'+str(y);
    return rets



def calchrhudu():
    sums=0
    a=0
    for i in listchrlength:
        sums=sums+int(listchrlength[a])
        a=a+1
    a=0
    dushu=0
    for i in listchrlength:
        dushu=dushu+(360/sums)*int(listchrlength[a])
        if 0<dushu<=90:
            if a==0:
                listhudu.append(dushu-2)
                listhudu.append(90-1)
                listhudu.append(90-dushu-1)
            else:
                listhudu.append((360/sums)*int(listchrlength[a])-2)
                listhudu.append(90-dushu+(360/sums)*int(listchrlength[a])-1)
                listhudu.append(90-dushu-1)
        if 90<dushu<=180:
            if 360-(dushu-90)+(360/sums)*int(listchrlength[a])>360:
                listhudu.append((360/sums)*int(listchrlength[a])-2)
                listhudu.append(360-(dushu-90)+(360/sums)*int(listchrlength[a])-360-1)
                listhudu.append(360-(dushu-90)-1)
            else:
                listhudu.append((360/sums)*int(listchrlength[a])-2)
                listhudu.append(360-(dushu-90)+(360/sums)*int(listchrlength[a])-1)
                listhudu.append(360-(dushu-90)-1)
        if 180<dushu<=270:
            listhudu.append((360/sums)*int(listchrlength[a])-2)
            listhudu.append(360-(dushu-90)+(360/sums)*int(listchrlength[a])-1)
            listhudu.append(360-(dushu-90)-1)
        if 270<dushu<361:
            listhudu.append((360/sums)*int(listchrlength[a])-2)
            listhudu.append(360-(dushu-90)+(360/sums)*int(listchrlength[a])-1)
            listhudu.append(360-(dushu-90)-1)
        a=a+1


    
def openlinkfile(filename):
    file2=open(filename,'r')
    contentf2=file2.read() 
    contsf2=contentf2.split('\n')
    a=0
    for i in contsf2:
        if len(contsf2[a])==0:
            break
        xinxi=contsf2[a].split('\t')
        chrname1=xinxi[0]
        namepos=listchrname.index(chrname1)
        pos1=xinxi[1]
        lenpos=listchrlength[namepos]
        jiaodu=listhudu[namepos*3+1]-listhudu[namepos*3]*int(pos1)/int(lenpos)
        if jiaodu<0:
            jiaodu=jiaodu+360
        calposbyangle(jiaodu,1)
        chrname2=xinxi[3]
        namepos=listchrname.index(chrname2)
    
        pos2=xinxi[4]
        lenpos=listchrlength[namepos]
        jiaodu=listhudu[namepos*3+1]-listhudu[namepos*3]*int(pos2)/int(lenpos)
        if jiaodu<0:
            jiaodu=jiaodu+360
        calposbyangle(jiaodu,1)
        a=a+1
    file2.close()


def opengenefile(filename):
    filegene=open(filename,'r')
    contentf2=filegene.read() 
    contsf2=contentf2.split('\n')
    a=0
    for i in contsf2:
        if len(contsf2[a])==0:
            break
        xinxi=contsf2[a].split('\t')
        chrname1=xinxi[0]
        namepos=listchrname.index(chrname1)
        pos1=xinxi[1]
        lenpos=listchrlength[namepos]
        jiaodu=listhudu[namepos*3+1]-listhudu[namepos*3]*int(pos1)/int(lenpos)
        if jiaodu<0:
            jiaodu=jiaodu+360
        calposbyangle(jiaodu,2)
        chrname2=xinxi[3]
        namepos=listchrname.index(chrname2)
    
        pos2=xinxi[4]
        lenpos=listchrlength[namepos]
        jiaodu=listhudu[namepos*3+1]-listhudu[namepos*3]*int(pos2)/int(lenpos)
        if jiaodu<0:
            jiaodu=jiaodu+360
        calposbyangle(jiaodu,2)
        a=a+1
    filegene.close()


def paintlinearcurve(colo):
    Path = mpath.Path
    a=0
    for i in range(0,int(len(list1)/4)):
        name='pp'+str(i+1)
        name = mpatches.PathPatch(Path([(list1[i+a], list1[i+1+a]), (0, 0), (list1[i+2+a], list1[i+3+a])],
             [Path.MOVETO, Path.CURVE3, Path.CURVE3]),fc="none", transform=ax.transData,color=colo,linewidth=0.5)
        ax.add_patch(name)
    
        a=a+3

def paintcolumn(yiwei,colour):
    a=0
    for i in listgene:
        if a+2>len(listgene):
            break
        angle = calc_angle(0,0,listgene[a],listgene[a+1])
    
        if listgene[a]>=0 and listgene[a+1]>0:
            y=yiwei*math.sin(math.radians(90-angle))
            x=(yiwei*yiwei-y*y)**0.5
            line1 = [(listgene[a]+x, listgene[a+1]+y), (listgene[a]+x+x*0.2, listgene[a+1]+y+y*0.2)]
            (line1_xs, line1_ys) = zip(*line1)
            ax.add_line(Line2D(line1_xs, line1_ys, linewidth=0.5, color=colour))

        elif listgene[a]<0 and listgene[a+1]<0:
            y=yiwei*math.sin(math.radians(90-angle))
            x=(yiwei*yiwei-y*y)**0.5
            line1 = [(listgene[a]-x, listgene[a+1]-y), (listgene[a]-x-x*0.2, listgene[a+1]-y-y*0.2)]
            (line1_xs, line1_ys) = zip(*line1)
            ax.add_line(Line2D(line1_xs, line1_ys, linewidth=0.5, color=colour))
        else:
            y=yiwei*math.sin(math.radians(angle))
            x=(yiwei*yiwei-y*y)**0.5
            if listgene[a]<0:
                line1 = [(listgene[a]-x, listgene[a+1]-y), (listgene[a]-x-x*0.2, listgene[a+1]-y-y*0.2)]
                (line1_xs, line1_ys) = zip(*line1)
                ax.add_line(Line2D(line1_xs, line1_ys, linewidth=0.5, color=colour))
            else:
                line1 = [(listgene[a]+x, listgene[a+1]+y), (listgene[a]+x+x*0.2, listgene[a+1]+y+y*0.2)]
                (line1_xs, line1_ys) = zip(*line1)
                ax.add_line(Line2D(line1_xs, line1_ys, linewidth=0.5, color=colour))
        a=a+2


def paintchr(cll):
    global outer_colors
    cmap = plt.get_cmap(cll)
    outer_colors = cmap(np.array([1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]))
    inner_colors = cmap(np.array([1, 2, 5, 6, 9, 10]))
    a=0
    yiweichr=widt*17.5
    for i in listchrname:
    
        lenpos=listchrlength[a]
        middle=listhudu[a*3+1]-listhudu[a*3]/2
        if middle<0:
            middle=middle+360
        rets=paintchrbyangle(middle+2,listchrname[a])
        xx=float(rets.split('\t')[0])
        yy=float(rets.split('\t')[1])
        angle = calc_angle(0,0,xx,yy)
    
        if xx>=0 and yy>0:
            angle = angle +5
            y=yiweichr*math.sin(math.radians(90-angle))
            x=(yiweichr*yiweichr-y*y)**0.5
            plt.text(xx+x, yy+y,listchrname[a],rotation=-angle+5,verticalalignment='center',horizontalalignment='center',color=outer_colors[a],size=7)
        elif xx<0 and yy<0:
            angle = angle +5
            y=yiweichr*math.sin(math.radians(90-angle))
            x=(yiweichr*yiweichr-y*y)**0.5
            plt.text(xx-x, yy-y,listchrname[a],rotation=-angle+5,verticalalignment='center',horizontalalignment='center',color=outer_colors[a],size=7)
        else:
            angle = angle -5
            y=yiweichr*math.sin(math.radians(angle))
            x=(yiweichr*yiweichr-y*y)**0.5
            if xx<0:
                ax.annotate(listchrname[a], xy=(xx-x, yy-y),rotation=angle-270+5,verticalalignment='center',horizontalalignment='center',color=outer_colors[a],size=7)
            else:
                ax.annotate(listchrname[a], xy=(xx+x, yy+y),rotation=angle-270+5,verticalalignment='center',horizontalalignment='center',color=outer_colors[a],size=7)
        a=a+1

    
def opengenedensityfile(filename):
    filegenedensity=open(filename,'r')
    contsgenedensity=filegenedensity.read().split('\n')
    a=0
    maxgenenum=0
    for i in contsgenedensity:
        if len(contsgenedensity[a])==0:
            break
        contssgenedensity=contsgenedensity[a].split('\t')
        if int(contssgenedensity[2])>int(maxgenenum):
            maxgenenum=contssgenedensity[2]
        a=a+1
    a=0
    for i in contsgenedensity:
        if len(contsgenedensity[a])==0:
            break
        xinxi=contsgenedensity[a].split('\t')
        chrname1=xinxi[0]
        namepos=listchrname.index(chrname1)
        pos1=xinxi[1]
        lenpos=listchrlength[namepos]
        jiaodu=listhudu[namepos*3+1]-listhudu[namepos*3]*int(pos1)/int(lenpos)
        if jiaodu<0:
            jiaodu=jiaodu+360
        calposbyangle(jiaodu,3)
        listgenedensity.append((int(xinxi[2])/(int(maxgenenum)+int(maxgenenum)*0.2))*0.2)
        a=a+1
    filegenedensity.close()


def paintgenedensityline(yiwei,colour):
    a=0
    for i in listgenedensity:
        if a+2>len(listgenedensity):
            break
        angle = calc_angle(0,0,listgenedensity[a],listgenedensity[a+1])
        if listgenedensity[a]>=0 and listgenedensity[a+1]>0:
            y=yiwei*math.sin(math.radians(90-angle))
            x=(yiwei*yiwei-y*y)**0.5
            line1 = [(listgenedensity[a]+x, listgenedensity[a+1]+y), (listgenedensity[a]+x+x*listgenedensity[a+2], listgenedensity[a+1]+y+y*listgenedensity[a+2])]
            (line1_xs, line1_ys) = zip(*line1)
            ax.add_line(Line2D(line1_xs, line1_ys, linewidth=0.5, color=colour))
        

        elif listgenedensity[a]<0 and listgenedensity[a+1]<0:
            y=yiwei*math.sin(math.radians(90-angle))
            x=(yiwei*yiwei-y*y)**0.5
            line1 = [(listgenedensity[a]-x, listgenedensity[a+1]-y), (listgenedensity[a]-x-x*listgenedensity[a+2], listgenedensity[a+1]-y-y*listgenedensity[a+2])]
            (line1_xs, line1_ys) = zip(*line1)
            ax.add_line(Line2D(line1_xs, line1_ys, linewidth=0.5, color=colour))
        
        else:
            y=yiwei*math.sin(math.radians(angle))
            x=(yiwei*yiwei-y*y)**0.5
            if listgenedensity[a]<0:
                line1 = [(listgenedensity[a]-x, listgenedensity[a+1]-y), (listgenedensity[a]-x-x*listgenedensity[a+2], listgenedensity[a+1]-y-y*listgenedensity[a+2])]
                (line1_xs, line1_ys) = zip(*line1)
                ax.add_line(Line2D(line1_xs, line1_ys, linewidth=0.5, color=colour))
            else:
                line1 = [(listgenedensity[a]+x, listgenedensity[a+1]+y), (listgenedensity[a]+x+x*listgenedensity[a+2], listgenedensity[a+1]+y+y*listgenedensity[a+2])]
                (line1_xs, line1_ys) = zip(*line1)
                ax.add_line(Line2D(line1_xs, line1_ys, linewidth=0.5, color=colour))
            
        a=a+3



def paintgenedensityscatter(yiwei,colour):
    a=0
    for i in listgenedensity:
        if a+2>len(listgenedensity):
            break
        angle = calc_angle(0,0,listgenedensity[a],listgenedensity[a+1])
        if listgenedensity[a]>=0 and listgenedensity[a+1]>0:
            y=yiwei*math.sin(math.radians(90-angle))
            x=(yiwei*yiwei-y*y)**0.5
            line1 = [(listgenedensity[a]+x, listgenedensity[a+1]+y), (listgenedensity[a]+x+x*listgenedensity[a+2], listgenedensity[a+1]+y+y*listgenedensity[a+2])]
            (line1_xs, line1_ys) = zip(*line1)
            ax.scatter([listgenedensity[a]+x+x*listgenedensity[a+2]], [listgenedensity[a+1]+y+y*listgenedensity[a+2]], c=colour,s=0.6)
        

        elif listgenedensity[a]<0 and listgenedensity[a+1]<0:
            y=yiwei*math.sin(math.radians(90-angle))
            x=(yiwei*yiwei-y*y)**0.5
            line1 = [(listgenedensity[a]-x, listgenedensity[a+1]-y), (listgenedensity[a]-x-x*listgenedensity[a+2], listgenedensity[a+1]-y-y*listgenedensity[a+2])]
            (line1_xs, line1_ys) = zip(*line1)
            ax.scatter([listgenedensity[a]-x-x*listgenedensity[a+2]], [listgenedensity[a+1]-y-y*listgenedensity[a+2]], c=colour,s=0.6)
        
        else:
            y=yiwei*math.sin(math.radians(angle))
            x=(yiwei*yiwei-y*y)**0.5
            if listgenedensity[a]<0:
                line1 = [(listgenedensity[a]-x, listgenedensity[a+1]-y), (listgenedensity[a]-x-x*listgenedensity[a+2], listgenedensity[a+1]-y-y*listgenedensity[a+2])]
                (line1_xs, line1_ys) = zip(*line1)
                ax.scatter([listgenedensity[a]-x-x*listgenedensity[a+2]], [listgenedensity[a+1]-y-y*listgenedensity[a+2]], c=colour,s=0.6)
            else:
                line1 = [(listgenedensity[a]+x, listgenedensity[a+1]+y), (listgenedensity[a]+x+x*listgenedensity[a+2], listgenedensity[a+1]+y+y*listgenedensity[a+2])]
                (line1_xs, line1_ys) = zip(*line1)
                plt.scatter([listgenedensity[a]+x+x*listgenedensity[a+2]], [listgenedensity[a+1]+y+y*listgenedensity[a+2]], c=colour,s=0.6)
            
        a=a+3




if __name__=="__main__":
    data1=''
    data2=''
    for opt_name,opt_value in opts:
        if opt_name in ('-i'):
            data1=opt_value
        if opt_name in ('-I'):
            data2=opt_value
fileconfig=open(data1,'r')
conts=fileconfig.read().split('\n')

a=0
t=0
linkfilename=''
chrcolo=''
linkradi=0
curvecolor=''
for i in conts:
    if len(conts[a])==0:
        break
    if conts[a][0]!='#':
        contss=conts[a].split('\t')
        if contss[0]=='chr':
            openchrfile(contss[1])
            chrcolo=contss[2]
            calchrhudu()
        if contss[0]=='link':
            linkfilename=contss[1]
            linkradi=widt*int(contss[3])
            curvecolor=contss[2]
        if contss[0]=='number':
            opengenefile(contss[1])
            paintcolumn(widt*int(contss[3]),contss[2])
            t=1
        if contss[0]=='densityline':
            opengenedensityfile(contss[1])
            paintgenedensityline(widt*int(contss[3]),contss[2])
            t=1
        if contss[0]=='densitydot':
            opengenedensityfile(contss[1])
            paintgenedensityscatter(widt*int(contss[3]),contss[2])
            t=1
    a=a+1
if t==0:
    radi=widt*22
    openlinkfile(linkfilename)
if t==1:
    openlinkfile(linkfilename)

paintlinearcurve(curvecolor)
paintchr(chrcolo)



ax.pie(listchrlength, radius=widt*22, colors=outer_colors,
       wedgeprops=dict(width=widt, edgecolor='w'),counterclock=False,startangle=90,explode=listexplode)


fileconfig.close

ax.spines["bottom"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["top"].set_visible(False)

ax.set(aspect="equal")
ax.axis([-3.5,3.5,-3.5,3.5])

ax.set_xticks([])

ax.set_yticks([])

plt.savefig('circos.png',dpi=800,figsize=(20,10))

plt.savefig('circos.svg',dpi=800)
plt.show()


