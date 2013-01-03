/*
 * rundelaunay.cpp
 *
 *  Created on: Aug 23, 2012
 *      Author: Trevor Kuhlengel
 */

#include "rundelaunay.h"
namespace std{

/** Executes a command line call given in the parameter
 \param cmd The full command that needs to be called
 \return The output from the command to stdout

 Code courtesy of the users of stackoverflow (http://stackoverflow.com/a/478960)
 */
std::string exec(char * cmd) {
	FILE* pipe = popen(cmd, "r");
	if (!pipe)
		return "ERROR";
	char buffer[130];
	std::string result = "";
	while (!feof(pipe)) {
		if (fgets(buffer, 128, pipe) != NULL)
			result += buffer;
	}
	pclose(pipe);
	return result;
}
/** Constructs the command line call to qDelaunay from the names provided
 \param params

 */
std::string buildCommand(string fileinput, string params) {
	char command[100]="";

	sprintf(command, "qdelaunay %s <%s", params.data(), fileinput.data());
	/*	}
	 else{
	 throw 21;
	 }
	 */
	return exec(command);
}

std::string buildAndPrint(int argc, char** arg) {
	//cout<<argc<<endl;
	std::string result("");
	bool clean=false;
	if (argc >= 2) {
		std::string str(arg[1]); //The modified filename for alteration
		std::string str2(".pdb");
		std::string filename //The original filename
					, coords;//The coordinates and
		filename.assign(str);
		int pos = str.find(str2);
		unsigned int npos=string::npos;

		//If the pdb extension is found, it needs to be converted to
		//A qhull readable format. This python script handles that lazily.
		if (pos != npos) {
			bool clean=true;
			str.replace(pos, str2.length(), ".qhc");
			char command[200]="";
			//This executes the python command and outputs to /tmp/<filename>.qhc
			sprintf(command, "python pdb2coord.py %s -s -o /tmp/%s \"reload\"",
					filename.data(), str.data() );

			//Assign the stdout to coords (if "-s" option is used, this contains radii)
			coords.assign(exec(command));
			filename.assign( "/tmp/").append(str);
		}
		string fil(filename);
		result = buildCommand(fil, "Qt i");
		if (clean){
			char rmCmd[50];
			sprintf(rmCmd, "rm -f %s \"reload\"", fil.data());
			system(rmCmd);
		}
		cout << result << endl;
		return result;
	}
	else{
		cout<<"Usage: Delaunay <filenameToRead>.<ext>\n"
			<<"\t Where <ext> can be \".pdb\" or any qhull readable format."
			<<endl;
		exit(0);
	}
	return result;
}
/******************************************************************************
 * All of the following functions are for string manipulation of inputs from
 * a certain python script. Once the script is no longer needed for input,
 * these functions can be deprecated.
 */
std::string getNext(std::string buffer, int &pos, bool eol){
	int endpos, npos=string::npos;
	std::string tok;
	if (eol){
		endpos=buffer.find("\n", pos);
	}
	else{
		endpos=buffer.find(" ", pos);
	}
	if (endpos != npos){
		tok.assign(buffer.substr(pos,endpos));
	}
	else{
		std::string buff;
		std::cout<<"Error in getInt: At position "<<pos<<"could not find "
				<<"matching separator for "<< (eol ? "return " : "space ");
		std::cout<<"character. Error code 90 thrown."<<endl;
		throw 90;
	}
	//int result=atoi(tok.data());
	pos=endpos;
	return tok;

}

Alpha::Sphere** loadCoords(std::string data, int &count, bool radii=true, double defaultradius=1.2){

	int  dim
		,pos=0 //current position
		,lpos
		,npos=string::npos;//Last position
	dim=atoi(getNext(data, pos, true).data());
	count=atoi(getNext(data, pos, true).data());
	Alpha::Vertex *pts=new Alpha::Vertex[count];
	Alpha::Sphere  **result=new Alpha::Sphere *[count];
	double *coords, radius;
	char line[50]="";
	coords=(radii ? new double[dim+1] : new double[dim]);
	double *vertPos;
	for (int i=0; i<count; i++){
		radius=1.2; //Default value if no radius is found
		//Find the end of the next line and send it to be parsed
		lpos=pos;
		pos=data.find("\n", lpos);
		coords=parsePoints(data.substr(lpos, pos),dim,radii);

		//Set the Vertex object's Values to the results
		vertPos=pts[i].getPos();
		for (int j=0; j<dim;j++){
			vertPos[j]=coords[j];
		}
		//Set the radius variable
		if (radii){
			radius=coords[dim];
		}
		result[i]=new Alpha::Sphere(&pts[i], radius);
	}
	return result;

}
double* parsePoints(std::string pyline ///The output from the python script pdb2coord using -s option
		,int dim ///Number of dimensions to the data (not including radius of the sphere), or, the number of doubles per line
		,bool radii ///Indicates if there is a radius column after the coordinates
		//Return for the radius of the sphere.
		){
	int pos=0, lpos, npos=string::npos;
	if (radii){dim++;}
	double *result=new double[dim];

	for(int i=0; i<dim; i++){
		lpos=pos;
		pos=pyline.find(" ", lpos);
		if(pos!=npos){
			result[i]=atof(pyline.substr(lpos,pos).data());
		}
		else if(npos!=(pos=pyline.find("\n", lpos))){
			result[i]=atof(pyline.substr(lpos,pos).data());
		}
		else{
			printf("Error in rundelaunay.cpp parsedoubles: Expected at least %i\
					 spaces but did not find number %i\n Line content: %s",
					 radii ? dim+1:dim, i, pyline.data());
			throw 71;
		}
	}
	return result;

}
set_t* linkSurfaces(int **links, int count, int dim){
	int i;

	set_t top = alloc_set(6);
	
	for(i = 0; i<count; i++){
		if(i%10000 == 0){ printf("Adjacency Info processed %i out of %i\n", i, count); }
		top = put_set(top, i);
		set_t inSet = alloc_set(0);

		inSet = tgrid->querySphere(links[3*i+0], links[3*i+1], links[3*i+2], .00001, inSet);

		top = associate_set(top, i, inSet);
	}
	// Go through triangles and add new things to top.
	/***********************************************************************************************
	***********************************************************************************************
	***********************************************************************************************/
	///populate adjacency info
	for(i = 0; i<numberTriangles; i++){
		set_t mySet;
		//Indices of the triangle
		int a = triangles[3*i+0];
		int b = triangles[3*i+1];
		int c = triangles[3*i+2];

		//make A the toplevel
		mySet = (set_t) mapsto_set(top, a);
		mySet = put_set(mySet, a); mySet = put_set(mySet, b); mySet = put_set(mySet, c);
		top = associate_set(top, a, mySet);

		mySet = (set_t) mapsto_set(top, b);
		mySet = put_set(mySet, a); mySet = put_set(mySet, b); mySet = put_set(mySet, c);
		top = associate_set(top, b, mySet);

		mySet = (set_t) mapsto_set(top, c);
		mySet = put_set(mySet, a); mySet = put_set(mySet, b); mySet = put_set(mySet, c);
		top = associate_set(top, c, mySet);
	}

}

}//Namespace std
