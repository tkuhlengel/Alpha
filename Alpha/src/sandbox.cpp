using namespace std;
int fill(bool **npdata, bool **output, int x, int y, int &xmax, int &ymax){
    if ((( x>=xmax) || (y >= ymax)) || ((y<0) || (x < 0))){
        return 0;
    }
    //Base Case
    else if(npdata[y,x]){
    	return 0;
    }
    else{
    	output[y][x]=true;
		//Up
        fill(npdata, output, x , (y + 1) , xmax, ymax);
        //left
        fill(npdata, output, x - 1 , y , xmax, ymax );
        //down
        fill(npdata, output, x , y - 1 , xmax, ymax);
        //right
        fill(npdata, output, x + 1, y , xmax, ymax);
        //If adding 3d, just add 2 z alterations
    }

    return 0;
}
int x=0, y=0, xm=xmax, ym=ymax;
fill(npdat,output, x,y,xm,ym);

