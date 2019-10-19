import cv2 
import numpy as np
import progressbar
import os
from matplotlib import pyplot as plt
#os.chdir('E:/1730895/Work/Abaqus/CT/')

def captureVideo(videoFile,dstDir='./',start=1,inc=0,imgPers=1,end=0,suffix=''):
    """
    inc为frame 增量
    imgPers表示每秒提取照片数目
    """
    if not os.path.isdir(os.path.abspath(dstDir)):
        print(os.path.abspath(dstDir))
        print(dstDir,' is not a valid directory')
        return

    cap = cv2.VideoCapture(videoFile)

    FPS=cap.get(cv2.CAP_PROP_FPS)
    FrameCount=cap.get(cv2.CAP_PROP_FRAME_COUNT)

    if start>FrameCount:
        start=1
    if inc==0:
        inc=np.ceil(FPS/imgPers)
    if end==0 or end>FrameCount:
        end=FrameCount

    extractCnt=int((end-start)/inc+1)
    print('Heigh: ',cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    print('Width: ',cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    print('FPS      : ',FPS)
    print('frame 总数:',FrameCount)
    print('From %d to %d with increment %d ,total %d pictures'%(start,end,inc,extractCnt))

    widgets = ['Progress: ',progressbar.Percentage(), ' ', progressbar.Bar('#'),' ', progressbar.Timer(),
            ' ', progressbar.ETA(), ' ', progressbar.FileTransferSpeed()]
    pbar=progressbar.ProgressBar(widgets=widgets,maxval=extractCnt)
    pbar.start()

    for i,frameIndex in enumerate(range(start,end+1,inc)):
        cap.set(cv2.CAP_PROP_POS_FRAMES,frameIndex)
        ret,frame=cap.read()
        if ret==False:
            break
        im=cv2.flip(frame,-1)
        cv2.imencode('.jpg',im)[1].tofile('%s/%s_%d.jpg'%(dstDir,suffix,frameIndex))
        pbar.update(i)
    cap.release()

def findThreshold(imageName):
    im=cv2.imread(imageName)
    #im=cv2.flip(im,-1)
    #gray=cv2.cvtColor(im,cv2.COLOR_BGR2GRAY)
    #cv2.imshow('gray',gray)
    #cv2.waitKey(0)
    #cv2.destroyAllWindows()

    img=cv2.cvtColor(im,cv2.COLOR_BGR2GRAY)
    ret,th1 = cv2.threshold(img,127,255,cv2.THRESH_BINARY)
    th2 = cv2.adaptiveThreshold(img,255,cv2.ADAPTIVE_THRESH_MEAN_C,cv2.THRESH_BINARY,11,2)
    th3 = cv2.adaptiveThreshold(img,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,cv2.THRESH_BINARY,11,2)
    titles = ['Original Image', 'Global Thresholding (v = 127)',
                'Adaptive Mean Thresholding', 'Adaptive Gaussian Thresholding']
    images = [img, th1, th2, th3]
    for i in range(4):
        plt.subplot(2,2,i+1),plt.imshow(images[i],'gray')
        plt.title(titles[i])
        plt.xticks([]),plt.yticks([])
    plt.show(block = False)

    cv2.namedWindow("Threshold")
    h,w=img.shape
    newShape=tuple(int(x/2) for x in  [w,h])
    cv2.resizeWindow("Threshold",newShape)
    def on_trace_bar_changed(pos):
        print(pos)
        cv2.setTrackbarPos('Track Bar', 'Threshold', pos)
        ret,th1 = cv2.threshold(img,pos,255,cv2.THRESH_BINARY)
        th2=cv2.resize(th1,newShape)
        cv2.imshow('Threshold',th2)
        while True:
            k=cv2.waitKey(0)
            if k==27 or k==ord('q'):
                cv2.destroyWindow('Threhold')
                print('Quit')
                break
            elif k==ord('s'):
                fnWithoutExt=os.path.splitext(imageName)[0]
                fname=fnWithoutExt+'-Threshold-%d.jpg'%pos
                cv2.imencode('.jpg',frame)[1].tofile(fname)
                print('image with threhold write to : ',fname)
                break
    
    cv2.createTrackbar('Track Bar', 'Threshold', 0, 255, on_trace_bar_changed)
    on_trace_bar_changed(78)
    
def gif():
    pass    
if __name__=="__main__":
    #os.chdir('E:/1730895/实验/')
    #videoFile='./基体-CT-断裂韧性-2019-08-21-21-56.MOV'
    #captureVideo(videoFile,'./基体-CT-断裂韧性-2019-08-21-21-56',start=900,inc=100,end=6200)
    #captureVideo(videoFile,'./基体-CT-断裂韧性-2019-08-21-21-56',start=6200,inc=50,end=7900)
    #captureVideo(videoFile,'./基体-CT-断裂韧性-2019-08-21-21-56',start=7900,inc=200,end=17430)
    #captureVideo(videoFile,'./基体-CT-断裂韧性-2019-08-21-21-56',start=17430,inc=1,end=17457)

    os.chdir('E:/1730895/实验/')
    videoFile='./基体-CT-断裂韧性-2019-08-21-21-56.MOV'
    captureVideo(videoFile,'E:/CT',start=900,inc=1000,end=6200,suffix='CT-2')
    captureVideo(videoFile,'E:/CT',start=6200,inc=340,end=7900,suffix='CT-2')
    captureVideo(videoFile,'E:/CT',start=7900,inc=1000,end=17430,suffix='CT-2')
    captureVideo(videoFile,'E:/CT',start=17430,inc=3,end=17457,suffix='CT-2')

    #videoFile='./基体-CT-断裂韧性-2019-08-06-15-52.mp4'
    #captureVideo(videoFile,r'E:/1730895/实验/基体-CT-断裂韧性-2019-08-06-15-52',start=7950,end=7980,inc=1)

    #findThreshold('./MQ-CT-1/CT_9000.jpg')