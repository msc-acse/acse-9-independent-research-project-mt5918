#define _USE_MATH_DEFINES

#include <algorithm>
#include <cassert>
#include <cmath>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <iomanip>
#include <iterator>
#include <map>
#include <fstream>
#include <functional>

#include <regex>
#include <string>
#include <vector>

#include <gmsh.h>

typedef std::string str;
typedef std::pair<int, int> intintpair;
typedef std::vector<int> intvec;
typedef std::pair<int, intvec> intintvecpair;
typedef std::map<int, int> intintmap;
typedef std::multimap<int, int> intintmultimap;		//Maps Node Numbers to Element Numbers 
typedef std::map<int, str> intstrmap;				//Maps Ints to Strings (Material Names)
typedef std::vector<str> strvec;
typedef std::vector<double> doublevec;
typedef intintmultimap::iterator intintmultimapit;

bool is_number(const std::string& s)
{
	//if there are non-digit characters, then string is not a number
	return !s.empty() && std::find_if(s.begin(),
		s.end(), [](char c) { return !isdigit(c); }) == s.end();
}

class Variables
{

};

class Material : public Variables
{

private:

	int number;
	str name;
	strvec data;

	double viscousdamping;
	double loadamplitude;
	double amplitude;
	double amplitudeinc;
	double loadamplitudeinc;
	double forcevelamplitudex;
	double forcevelamplitudey;

	double density;
	double youngsmod;
	double nu;
	double massdamping;
	double elasticpenalty;
	double contactpenalty;
	double modeienergyrate;
	double modeiienergyrate;
	double tensilestrength;
	double internalfriction;
	double internalcohesion;
	double porefluidpressure;
	double jointfriction;
	double jrc0;
	double jcs0;
	double jointsamplesize;
	double interfacefriction;
	str problems2d;

	double lambda;
	double mu;
	int fractureflag;
	int numberofmeshrefinements;
	int propertytype;
	double slidingfriction;
	int jointproperty;
	double dtcrit;

public:
	Material() {};   // This is the constructor declaration
	~Material() {};  // This is the destructor: declaration

	int getnumber() { return this->number; }
	str getname() { return this->name; }
	strvec getdata() { return this->data; }

	double getviscousdamping() { return this->viscousdamping; }
	double getloadamplitude() { return this->loadamplitude; }
	double getamplitude() { return this->amplitude; }
	double getamplitudeinc() { return this->amplitudeinc; }
	double getloadamplitudeinc() { return this->loadamplitudeinc; }
	double getforcevelamplitudex() { return this->forcevelamplitudex; }
	double getforcevelamplitudey() { return this->forcevelamplitudey; }

	double getdensity() { return this->density; }
	double getyoungsmod() { return this->youngsmod; };
	double getnu() { return this->nu; };
	double getmassdamping() { return this->massdamping; };
	double getelasticpenalty() { return this->elasticpenalty; };
	double getcontactpenalty() { return this->contactpenalty; };
	double getmodeienergyrate() { return this->modeienergyrate; };
	double getmodeiienergyrate() { return this->modeiienergyrate; };
	double gettensilestrength() { return this->tensilestrength; };
	double getinternalfriction() { return this->internalfriction; };
	double getinternalcohesion() { return this->internalcohesion; };
	double getporefluidpressure() { return this->porefluidpressure; };
	double getjointfriction() { return this->jointfriction; };
	double getjrc0() { return this->jrc0; };
	double getjcs0() { return this->jcs0; };
	double getjointsamplesize() { return this->jointsamplesize; };
	double getinterfacefriction() { return this->interfacefriction; };
	str getproblems2d() { return this->problems2d; };

	double getlambda() { return this->lambda; }
	double getmu() { return this->mu; }
	int getfractureflag() { return this->fractureflag; }
	int getnumberofmeshrefinements() { return this->numberofmeshrefinements; }
	int getpropertytype() { return this->propertytype; }
	double getslidingfriction() { return this->slidingfriction; }
	int getjointproperty() { return this->jointproperty; }
	double getdtcrit() { return this->dtcrit; }

	void setnumber(int number) { this->number = number; }
	void setname(str name) { this->name = name; }
	void setdata(strvec data) { this->data = data; }

	void setviscousdamping(double viscousdamping) { this->viscousdamping = viscousdamping; }
	void setloadamplitude(double loadamplitude) { this->loadamplitude = loadamplitude; }
	void setamplitude(double amplitude) { this->amplitude = amplitude; }
	void setamplitudeinc(double amplitudeinc) { this->amplitudeinc = amplitudeinc; }
	void setloadamplitudeinc(double loadamplitude) { this->loadamplitude = loadamplitude; }

	void setdensity(double density) { this->density = density; }
	void setyoungsmod(double youngsmod) { this->youngsmod = youngsmod; }
	void setnu(double nu) { this->nu = nu; }
	void setmassdamping(double massdamping) { this->massdamping = massdamping; }
	void setelasticpenalty(double elasticpenalty) { this->elasticpenalty = elasticpenalty; }
	void setcontactpenalty(double contactpenalty) { this->contactpenalty = contactpenalty; }
	void setmodeienergyrate(double modeienergyrate) { this->modeienergyrate = modeienergyrate; }
	void setmodeiienergyrate(double modeiienergyrate) { this->modeiienergyrate = modeiienergyrate; }
	void settensilestrength(double tensilestrength) { this->tensilestrength = tensilestrength; }
	void setinternalfriction(double internalfriction) { this->internalfriction = internalfriction; }
	void setinternalcohesion(double internalcohesion) { this->internalcohesion = internalcohesion; }
	void setporefluidpressure(double porefluidpressure) { this->porefluidpressure = porefluidpressure; }
	void setjointfriction(double jointfriction) { this->jointfriction = jointfriction; }
	void setjrc0(double jrc0) { this->jrc0 = jrc0; }
	void setjcs0(double jcs0) { this->jcs0 = jcs0; }
	void setjointsamplesize(double jointsamplesize) { this->jointsamplesize = jointsamplesize; }
	void setinterfacefriction(double interfacefriction) { this->interfacefriction = interfacefriction; }
	void setproblems2d(str problems2d) { this->problems2d = problems2d; }

	void setlambda(double lambda) { this->lambda = lambda; }
	void setmu(double mu) { this->mu = mu; }
	void setfractureflag(int fractureflag) { this->fractureflag = fractureflag; }
	void setnumberofmeshrefinements(int numberofmeshrefinements) { this->numberofmeshrefinements = numberofmeshrefinements; }
	void setpropertytype(int propertytype) { this->propertytype = propertytype; }
	void setslidingfriction(double slidingfriction) { this->slidingfriction = slidingfriction; }
	void setjointproperty(int jointproperty) { this->jointproperty = jointproperty; }
	void setdtcrit(double dtcrit) { this->dtcrit = dtcrit; }

};

typedef std::map<str, Material> strmatmap;			//Maps Strings to Materials
typedef std::map<int, Material> intmatmap;			//Maps Ints to Materials

class PropertyVariables
{
	friend class Material;

private:
	int maxnumprop;
	int numprop;
	strmatmap strmat;
	intmatmap intmat;
	intstrmap matname;
	int nummats;
	doublevec dtcrits;
	double dtcrit;

	double minelev;
	double minedge;

	double projectilemeshsize;
	double plymeshsize;
	double interlayermeshsize;

public:
	PropertyVariables() {};   // This is the constructor declaration
	~PropertyVariables() {};  // This is the destructor: declaration

	int getmaxnumprop() { return this->maxnumprop; }
	int getnumprop() { return this->numprop; }
	Material getstrmat(str name);
	Material getintmat(int number);
	str getmatname(int number);
	int getmatnum(str name);

	int getnummats() { return this->nummats; }
	double getdtcrit();
	doublevec getdtcrits();

	double getminelev() { return this->minelev; }
	double getminedge() { return this->minedge; }

	double getprojectilemeshsize() { return this->projectilemeshsize; }
	double getplymeshsize() { return this->plymeshsize; };
	double getinterlayermeshsize() { return this->interlayermeshsize; }

	void setmaxnumprop(int maxnumprop) { this->maxnumprop = maxnumprop; }
	void setnumprop(int numprop) { this->numprop = numprop; }

	void setstrmat(str name, Material Mat) { this->strmat[name] = Mat; }
	void setintmat(int number, Material Mat) { this->intmat[number] = Mat; }
	void setmatname(int number, str name) { this->matname[number] = name; }

	void setnummats(int nummats) { this->nummats = nummats; }
	void setdtcrit(double dtcrit) { this->dtcrit = dtcrit; }
	void setdtcrits(doublevec dtcrits) { this->dtcrits = dtcrits; }

	void setminelev(double minelev) { this->minelev = minelev; }
	void setminedge(double minedge) { this->minedge = minedge; }

	void setprojectilemeshsize(double projectilemeshsize) { this->projectilemeshsize = projectilemeshsize; }
	void setplymeshsize(double plymeshsize) { this->plymeshsize = plymeshsize; }
	void setinterlayermeshsize(double interlayermeshsize) { this->interlayermeshsize = interlayermeshsize; }
};

int PropertyVariables::getmatnum(str name)
{
	//find the material id from the material name, e.g. "steel" ~> 1

	for (intstrmap::iterator it = this->matname.begin(); it != this->matname.end(); it++) //iterate through material intstrmap
	{
		if (it->second.compare(name)) //if corresponding string is found
		{
			return it->first; //return matid
		}
	}

	std::cerr << "Error (in getmatnum): Material not set; Corresponding material id not found. Check your mat file." << std::endl;
	return -1;
}

str PropertyVariables::getmatname(int number)
{
	//get name of material from material id, e.g. 1 ~> "steel"
	for (intstrmap::iterator it = this->matname.begin(); it != this->matname.end(); it++) //iterate through intstrmap
	{
		if (it->first == number) //if corresponding mat id is found
		{
			return it->second; //return mat name
		}
	}

	std::cerr << "Error (in getmatname): string index out of bound; material name not found. Check your mat file." << std::endl;
	return " ";
}

Material PropertyVariables::getintmat(int number)
{
	//get material object from material id, e.g. 1 ~> Material object for steel

	for (intmatmap::iterator it = this->intmat.begin(); it != this->intmat.end(); it++) //iterate through intmatmap
	{
		if (it->first == number) //if mat id corresponds
		{
			return it->second; //return the material object
		}
	}

	std::cerr << "Error: material id out of bound; input id: " << number << ". Check your mat file." << std::endl;
	Material Error;
	return Error;
}

Material PropertyVariables::getstrmat(str name)
{
	//get material object from material name, e.g. "steel" ~> material object for steel
	for (strmatmap::iterator it = this->strmat.begin(); it!= this->strmat.end(); it++)
	{
		if (it->first == name)
		{
			return it->second;
		}
	}

	std::cerr << "Error: Material index out of bound. Check your mat file." << std::endl;
	Material Error;
	return Error;
}

doublevec PropertyVariables::getdtcrits()
{
	//calculate the critical time steps from the materials;

	this->dtcrits.clear(); //clear the vector
	for (strmatmap::iterator it = this->strmat.begin(); it != this->strmat.end(); ++it) //iterate through the materials
	{ 
		assert(it->second.getdtcrit() > 0 && "Error (in getdtcrits): critical time step is invalid.\n");
		this->dtcrits.push_back(it->second.getdtcrit());  //store the critical time steps from the materials
	}
	return this->dtcrits;
}


double PropertyVariables::getdtcrit()
{
	//calculate the total (smallest) critical time step from the vector of critical time steps
	this->dtcrits = this->getdtcrits(); //only loop once
	this->dtcrit = this->dtcrits[std::distance(this->dtcrits.begin(), std::min_element(this->dtcrits.begin(), this->dtcrits.end()))];
	assert(this->dtcrit > 0 && "Error (in getdtcrit): critical time step is invalid.\n");

	return this->dtcrit;
}

class ProblemVariables: public Variables
{

private:

	strvec problemdata;

	double realsimtime;
	int numberoutfiles;
	double maxndt;
	double gravity;
	double dtcrit;
	double dt;
	double outfreq;
	double restartfreq;
	int gravitysettlingstage;
	int loadrampingstage;
	int maxdim;
	double maxforce;
	int maxv;
	double maxstress;
	int maxdispl;
	double maxjointaperture;
	int maxcontcouples;
	double buffersize;
	str accuracy;
	str fricmodel;
	str initialaperturecorr;
	str genoutonly;

	int currentndt;
	double coordsize;
	double currentt;
	int contcouples;
	int numchars;
	int currentit;

	double simtime;
	int numoutfiles;

public:
	ProblemVariables() 
	{
		this->setcurrentit(0);
		this->setnumchars(4);
	};   // This is the constructor declaration
	~ProblemVariables() {};  // This is the destructor: declaration

	strvec getproblemdata() { return this->problemdata; }

	str getaccuracy() { return this->accuracy; }
	str getfricmodel() { return this->fricmodel; }
	str getinitialaperturecorr() { return this->initialaperturecorr; }

	double getrealsimtime() { return this->realsimtime; }
	int getnumberoutfiles() { return this->numberoutfiles; }
	double getmaxndt() { return this->maxndt; }
	double getgravity() { return this->gravity; }
	double getdt() { return this->dt; }
	double getdtcrit() { return this->dtcrit; }
	double getoutfreq() { return this->outfreq; }
	double getrestartfreq() { return this->restartfreq; }
	int getgravitysettlingstage() { return this->gravitysettlingstage; }
	int getloadrampingstage() { return this->loadrampingstage; }
	int getmaxdim() { return this->maxdim; }
	double getmaxforce() { return this->maxforce; }
	int getmaxv() { return this->maxv; }
	double getmaxstress() { return this->maxstress; }
	int getmaxdispl() { return this->maxdispl; }
	double getmaxjointaperture() { return this->maxjointaperture; }
	int getmaxcontcouples() { return this->maxcontcouples; }
	double getbuffersize() { return this->buffersize; }

	double getcurrentdt() { return this->currentndt; }
	double getcoordsize() { return this->coordsize; }
	double getcurrentt() { return this->currentt; }
	int getcontcouples() { return this->contcouples; }

	str getgenoutonly() { return this->genoutonly; }
	int getcurrentit() { return this->currentit; }
	int getnumchars() { return this->numchars; }

	void setrealsimtime(double realsimtime) { this->realsimtime = realsimtime; }
	void setproblemdata(strvec problemdata) { this->problemdata = problemdata; }

	void setnumberoutfiles(int numberoutfiles) { this->numberoutfiles = numberoutfiles; }
	void setmaxndt(double maxndt) { this->maxndt = maxndt; }
	void setgravity(double gravity) { this->gravity = gravity; }
	void setdt(double dt) { this->dt = dt; }
	void setdtcrit(double dtcrit) { this->dtcrit = dtcrit; }
	void setoutfreq(double outfreq) { this->outfreq = outfreq; }
	void setrestartfreq(double restartfreq) { this->restartfreq = restartfreq; }
	void setgravitysettlingstage(int gravitysettlingstage) { this->gravitysettlingstage = gravitysettlingstage; }
	void setloadrampingstage(int loadrampingstage) { this->loadrampingstage = loadrampingstage; }
	void setmaxdim(int maxdim) { this->maxdim = maxdim; }
	void setmaxforce(double maxforce) { this->maxforce = maxforce; }
	void setmaxv(int maxv) { this->maxv = maxv; }
	void setmaxstress(double maxstress) { this->maxstress = maxstress; }
	void setmaxdispl(int maxdispl) { this->maxdispl = maxdispl; }
	void setmaxjointaperture(double maxjointaperture) { this->maxjointaperture = maxjointaperture; }
	void setmaxcontcouples(int maxcontcouples) { this->maxcontcouples = maxcontcouples; }
	void setbuffersize(double buffersize) { this->buffersize = buffersize; }
	void setaccuracy(str accuracy) { this->accuracy = accuracy; }
	void setfricmodel(str fricmodel) { this->fricmodel = fricmodel; }
	void setinitialaperturecorr(str initialaperturecorr) { this->initialaperturecorr = initialaperturecorr; }

	void setcurrentndt(int currentndt) { this->currentndt = currentndt; }
	void setcoordsize(int coordsize) { this->coordsize = coordsize; }
	void setcurrentt(double currentt) { this->currentt = currentt; }
	void setcontcouples(int contcouples) { this->contcouples = contcouples; }

	void setgenoutonly(str genoutonly) { this->genoutonly = genoutonly; }
	void setcurrentit(int currentit) { this->currentit = currentit; }
	void setnumchars(int numchars) { this->numchars = numchars; }

	double getsimtime() { return this->simtime; }
	int getnumoutfiles() { return this->numoutfiles; }

	void setsimtime(double simtime) { this->simtime = simtime; }
	void setnumoutfiles(int numoutfiles) { this->numoutfiles = numoutfiles; }
};


class UserVariables : public Variables
{
private:

	strvec userdata;

	str infile;
	str makenewfile;
	str outfile;
	str login;
	str hpcdir;
	str localdir;
	str qsubdir;
	str windowsmayavidir;
	str bashmayavidir;
	str runlocal;
	str windowsgiddir;
	str bashgiddir;
	str windowsydir;
	str bashydir;
	str windowsparaviewdir;
	str bashparaviewdir;

public:
	UserVariables() {};   // This is the constructor declaration
	~UserVariables() {};  // This is the destructor: declaration

	strvec getuserdata() { return this->userdata; }
	str getmakenewfile() { return this->makenewfile; }
	str getoutfile() { return this->outfile; }
	str getinfile() { return this->infile; }
	str getlogin() { return this->login; }
	str gethpcdir() { return this->hpcdir; }
	str getlocaldir() { return this->localdir; }
	str getqsubdir() { return this->qsubdir; }
	str getwindowsmayavidir() { return this->windowsmayavidir; }
	str getbashmayavidir() { return this->bashmayavidir; }
	str getrunlocal() { return this->runlocal; }
	str getwindowsgiddir() { return this->windowsgiddir; }
	str getbashgiddir() { return this->bashgiddir; }
	str getwindowsydir() { return this->windowsydir; }
	str getbashydir() { return this->bashydir; }
	str getwindowsparaviewdir() { return this->windowsparaviewdir; }
	str getbashparaviewdir() { return this->bashparaviewdir; }

	void setuserdata(strvec userdata) { this->userdata = userdata; }

	void setmakenewfile(str makenewfile) { this->makenewfile = makenewfile; }
	void setoutfile(str outfile) { this->outfile = outfile; }
	void setinfile(str infile) { this->infile = infile; }
	void setlogin(str login) { this->login = login; }
	void sethpcdir(str hpcdir) { this->hpcdir = hpcdir; }
	void setlocaldir(str localdir) { this->localdir = localdir; }
	void setqsubdir(str qsubdir) { this->qsubdir = qsubdir; }
	void setwindowsmayavidir (str windowsmayavidir) { this->windowsmayavidir = windowsmayavidir; }
	void setbashmayavidir(str  bashmayavidir) { this->bashmayavidir = bashmayavidir; }
	void setrunlocal(str  runlocal) { this->runlocal = runlocal; }
	void setwindowsgiddir(str  windowsgiddir) { this->windowsgiddir = windowsgiddir; }
	void setbashgiddir(str  bashgiddir) { this->bashgiddir = bashgiddir; }
	void setwindowsydir(str  windowsydir) { this->windowsydir = windowsydir; }
	void setbashydir(str  bashydir) { this->bashydir = bashydir; }
	void setwindowsparaviewdir(str  windowsparaviewdir) { this->windowsparaviewdir = windowsparaviewdir; }
	void setbashparaviewdir(str  bashparaviewdir) { this->bashparaviewdir = bashparaviewdir; }
};

class Node
{
private:
	int number;
	double dx;
	double dy;
	double vx;
	double vy;
	double fx;
	double fy;
	int matid;
	double fcx; //force constraint in x
	double fcy;
	double vcx; //velocity constraint
	double vcy;
	double dcx; //displacement constraint
	double dcy;
	int cflag; //constraintflag
	int initflag; //initial data flag;
	int boundaryflag;

public:
	int getnumber() { return this->number; }

	double getdx() { return this->dx; }
	double getdy() { return this->dy; }

	double getvx() { return this->vx; }
	double getvy() { return this->vy; }

	double getfx() { return this->fx; }
	double getfy() { return this->fy; }

	int getmatid() { return this->matid; }
	int getcflag() { return this->cflag; }

	double getdcx() { return this->dcx; }
	double getdcy() { return this->dcy; }

	double getfcx() { return this->fcx; }
	double getfcy() { return this->fcy; }

	double getvcx() { return this->vcx; }
	double getvcy() { return this->vcy; }

	int getinitflag() { return this->initflag; }
	int getboundaryflag() { return this->boundaryflag; }
	///////////////////////////////////
	///////////////////////////////////

	void setnumber(int number) { this->number = number; }

	void setdx(double dx) { this->dx = dx; }
	void setdy(double dy) { this->dy = dy; }

	void setvx(double vx) { this->vx = vx; }
	void setvy(double vy) { this->vy = vy; }

	void setfx(double fx) { this->fx = fx; }
	void setfy(double fy) { this->fy = fy; }

	void setmatid(int matid) { this->matid = matid; }
	void setcflag(int cflag) { this->cflag = cflag; }

	void setdcx(double dcx) { this->dcx = dcx; }
	void setdcy(double dcy) { this->dcy = dcy; }

	void setfcx(double fcx) { this->fcx = fcx; }
	void setfcy(double fcy) { this->fcy = fcy; }

	void setvcx(double vcx) { this->vcx = vcx; }
	void setvcy(double vcy) { this->vcy = vcy; }

	void setinitflag(int flag) { this->initflag = flag; }
	void setboundaryflag(int flag) { this->boundaryflag = flag; }
};

class Element
{
private:
	int number;
	int fractureflag;
	int fracturetype; //Fracture or Visco-Elastic
	intintmap nodenums; //node numbers of elements
	int numnodes; // number of nodes in element

public:

	int getnumber() { return this->number; }
	int getfractureflag() { return this->fractureflag; }
	int getfracturetype() { return this->fracturetype; }
	int getnodenum(int counter) { return this->nodenums[counter]; }
	int getnumnodes() { return this->numnodes; }

	void setnumber(int elementnumber) { this->number = elementnumber; }
	void setfractureflag() { this->fractureflag = fractureflag; }
	void setfracturetype() { this->fracturetype = fracturetype; }
	void setnodenum(int counter,int nodenumber) { this->nodenums[counter] = nodenumber; }
	void setnumnodes(int numnodes) { this->numnodes = numnodes; }
};

typedef std::map<int, Node> intnodemap;

class Surface
{
private:
	int matid; //e.g. 0 (steel)
	int number; // e.g. 0, 1, 2, 3, 4
	str type; //e.g. circle, rectangle
	str matname;
	intnodemap intnodes;

public:
	int getmatid() { return this->matid; }
	int getnumber() { return this->number; }
	str gettype() { return this->type; }
	str getmatname() { return this->matname; }
	Node getintnode(int nodenumber) { return this->intnodes[nodenumber]; }

	void setmatid(int matid) { this->matid = matid; }
	void setnumber(int number) { this->number = number; }
	void settype(str type) { this->type = type; }
	void setmatname(str matname) { this->matname = matname; }
	void setintnode(int nodenumber, Node N) { this->intnodes[nodenumber] = N; }
};


typedef std::map<int, Element> intelemap;
typedef std::pair<int, Element> intelepair;
typedef std::pair<int, Node> intnodepair;
typedef std::map<int, Surface> intsurfacemap;

class GeometryVariables : public Variables
{

private:

	int dim;
	int maxdim;
	int maxnumele;
	int numele;
	int maxstatevar;
	int statevar;
	int maxnodesperfe;
	int nodesperfe;
	int numnodes;
	int maxnumnodes;
	int numconstraints;
	int numinit; //number of initial conditions
	int numboundary;
	int numsurfaces;

	double minelev; //minimum element volume
	double minedge; //minimum edge

	strvec elementdata;
	strvec nodedata;
	strvec initnodes;
	strvec initdata;
	strvec nforce;
	strvec nvel;

	intintmultimap nodeelemap; //maps nodes to elements
	intelemap intelemap;
	intnodemap intnodemap;
	intsurfacemap intsurfacemap;

public:
	GeometryVariables() 
	{
		this->dim = -1;
		this->maxdim = -1;
		this->maxnumele = -1;
		this->numele = -1;
		this->maxstatevar = -1;
		this->statevar = -1;
		this->maxnodesperfe = -1;
		this->nodesperfe = -1;
		this->numnodes = -1;
		this->maxnumnodes = -1;
		this->numboundary = -1;

		this->numinit = 0;//number of initial conditions
		this->numconstraints = 0;
	};

	~GeometryVariables() {};

	int getdim() { return this->dim; }
	int getmaxdim() { return this->maxdim; }

	double getminelev() { return this->minelev; }
	double getminedge() { return this->minedge; }

	strvec getelementdata() { return this->elementdata; }
	strvec getnodes() { return this->nodedata; }
	strvec getinitnodedata() { return this->initnodes; }
	strvec getnforce() { return this->nforce; }
	strvec getnvel() { return this->nvel; }

	int getmaxnumele() { return this->maxnumele; }
	int getnumele() { return this->numele; }
	int getmaxstatevar() { return this->maxstatevar; }
	int getstatevar() { return this->statevar; }
	int getmaxnodesperfe() { return this->maxnodesperfe; }
	int getnodesperfe() { return this->nodesperfe; }
	int getnumnodes() { return this->numnodes; }
	int getmaxnumnodes() { return this->maxnumnodes; }
	int getnumconstraints() { return this->numconstraints; }
	int getnuminit() { return this->numinit; }
	int getnumboundary() { return this->numboundary; }
	int getnumsurfaces() { return this->numsurfaces; }

	intvec getnodeele(int nodenumber);
	Element getintele(int elenumber);
	Node getintnode(int nodenumber);
	Surface getintsurface(int surfacenumber);

	void setdim(int dim) { this->dim = dim; }
	void setmaxdim(int maxdim) { this->maxdim = maxdim; }

	void setminelev(double minelev) { this->minelev = minelev; }
	void setminedge(double minedge) { this->minedge = minedge; }

	void setmaxnumele(int maxnumele) { this->maxnumele = maxnumele; }
	void setnumele(int numele) { this->numele = numele; }
	void setmaxstatevar(int maxstatevar) { this->maxstatevar = maxstatevar; }
	void setstatevar(int statevar) { this->statevar = statevar; }
	void setmaxnodesperfe(int maxnodesperfe) { this->maxnodesperfe = maxnodesperfe; }
	void setnodesperfe(int nodesperfe) { this->nodesperfe = nodesperfe; }
	void setnumnodes(int numnodes) { this->numnodes = numnodes; }
	void setmaxnumnodes(int maxnumnodes) { this->maxnumnodes = maxnumnodes; }
	void setnuminit(int setnuminit) { this->numinit = numinit; }

	void setnumconstraints(int numconstraints) { this->numconstraints = numconstraints; }
	void setelementdata(strvec elements) { this->elementdata = elements; }
	void setnodedata(strvec nodes) { this->nodedata = nodes; }
	void setinitnodes(strvec initnodes) { this->initnodes = initnodes; }
	void setnforce(strvec nforce) { this->nforce = nforce; }
	void setnvel(strvec nvel) { this->nvel = nvel; }

	void setnodeele(int nodenum, int elenum);
	void setintnode(int nodenumber, Node N) { this->intnodemap[nodenumber] = N; }
	void setintele(int elenumber, Element E) { this->intelemap[elenumber] = E; }
	void setnumboundary(int numboundary) { this->numboundary = numboundary; }

	void setintsurface(int surfacenum, Surface S) { this->intsurfacemap[surfacenum] = S; }
	void setnumsurfaces(int numsurfaces) { this->numsurfaces = numsurfaces; }
};

Surface GeometryVariables::getintsurface(int surfacenumber)
{
	return this->intsurfacemap[surfacenumber];
}

Element GeometryVariables::getintele(int elenumber)
{
	assert(0 <= elenumber <= this->getnumele() && "Error (in getintele): Element number out of bound.\n");
	return this->intelemap[elenumber];
}

Node GeometryVariables::getintnode(int nodenumber)
{
	assert(0 <= nodenumber <= this->getnumnodes() && "Error (in getintnode): Node number out of bound.\n");
	return this->intnodemap[nodenumber];
}

void GeometryVariables::setnodeele(int nodenum, int elenum)
{
	assert(0 <= nodenum < this->getmaxnumnodes() && "Error(in setnodeele) : Node number invalid.\n");
	assert(0 <= elenum < this->getmaxnumele() && "Error(in setnodeele) : Element number invalid.\n");

	std::pair<intintmultimapit, intintmultimapit> elements = nodeelemap.equal_range(nodenum);
	int count = std::distance(elements.first, elements.second);

	assert(count < this->getmaxnodesperfe() && "Error (in setnodeele): Maximum number of nodes per FE exceeded.\n");

	this->nodeelemap.insert(intintpair (nodenum, elenum));
}

intvec GeometryVariables::getnodeele(int nodenum)
{
	assert(0 <= nodenum < this->getmaxnumnodes() && "Error(in getnodeele) : Node number invalid.\n");

	intvec elenums;
	std::pair<intintmultimapit, intintmultimapit> elements = nodeelemap.equal_range(nodenum);
	for (intintmultimap::iterator it = elements.first; it != elements.second; it++)
	{
		elenums.push_back(it->second);
	}
	return elenums;
}

//intvec GeometryVariables::getnodenums(int elenum)
//{
//	std::pair<intintmultimapit, intintmultimapit> nodes = nodeelemap.equal_range(elenum);
//	for (intintmultimap::iterator it = elements.first; it != elements.second; it++)
//		std::cout << it->second << std::endl;
//}

class Customise
{
	friend GeometryVariables;
	friend Node;
private:
	str projectilemat;
	strvec data;
	str userfile;
	str matfile;
	str ctrlfile;
	str nodefile;
	str nodevfile;
	str nodeffile;
	str ndcfile;
	str nvcfile;
	str nfcfile;
	str elementfile;
	str initfile;
	str geometryfile;

public:

	str getuserfile() { return this->userfile; }
	str getmatfile() { return this->matfile; }
	str getctrlfile() { return this->ctrlfile; }
	str getgeometryfile() { return this->geometryfile; }

	str getnodefile() { return this->nodefile; }
	str getnodevfile() { return this->nodevfile; }
	str getnodeffile() { return this->nodeffile; }
	str getndcfile() { return this->ndcfile; }
	str getnvcfile() { return this->nvcfile; }
	str getnfcfile() { return this->nfcfile; }
	str getelementfile() { return this->elementfile; }
	str getinitfile() { return this->initfile; }

	strvec getdata() { return this->data; }

	//void setprojectilenodesmat(str name, PropertyVariables &Property, GeometryVariables &Geometry);
	void setdata(strvec data) { this->data = data; }

	void setuserfile(str userfile) { this->userfile = userfile;}
	void setmatfile(str matfile) { this->matfile = matfile; }
	void setctrlfile(str ctrlfile) { this->ctrlfile = ctrlfile; }
	void setnodefile(str nodefile) { this->nodefile = nodefile; }
	void setnodevfile(str nodevfile) { this->nodevfile = nodevfile; }
	void setnodeffile(str nodeffile) { this->nodeffile = nodeffile; }
	void setndcfile(str ndcfile) { this->ndcfile = ndcfile; }
	void setnvcfile(str nvcfile) { this->nvcfile = nvcfile; }
	void setnfcfile(str nfcfile) { this->nfcfile = nfcfile; }
	void setelementfile(str elementfile) { this->elementfile = elementfile; }
	void setinitfile(str initfile) { this->initfile = initfile; }
	void setgeometryfile(str geometryfile) { this->geometryfile = geometryfile; }
};

//void Customise::setprojectilenodesmat(str name, PropertyVariables &Property, GeometryVariables &Geometry)
//{
//	this->projectilemat = name;
//	for (int i = 0; i != Geometry.getnumnodes(); i++)
//	{
//		if (Geometry.getintnode(i).getvx() != 0) //proectile material defined by nodes moving in x direction
//		{
//			std::cout << "Property.getmatnum(this->projectilemat): " << Property.getmatnum(this->projectilemat) << std::endl;
//			Geometry.getintnode(i).setmatid(Property.getmatnum(this->projectilemat));
//		}
//	}
//}

class Input
{
public:
	void readgeometry(GeometryVariables &Geometry, Surface &S, Node &N, Customise &Custom);
	void readelements(GeometryVariables &Geometry, Element &E, Customise &Custom);
	void readnodes(GeometryVariables &Geometry, Node &N, Customise &Custom, PropertyVariables &Property);
	void readuser(UserVariables &User, Customise &Custom);
	void readproblem(ProblemVariables &Problem, PropertyVariables &Property, Customise &Custom, GeometryVariables &Geometry);
	void readmaterial(Material Mat, PropertyVariables &Property, ProblemVariables &Problem, Customise &Custom, GeometryVariables &Geometry);
	void readnodalvelocities(GeometryVariables &Geometry, Node &N, Customise &Custom);
	void readnodalforces(GeometryVariables &Geometry, Node &N, Customise &Custom);
	void readdisplconstraints(GeometryVariables &Geometry, Node &N, Customise &Custom);
	void readforceconstraints(GeometryVariables &Geometry, Node &N, Customise &Custom);
	void readvelconstraints(GeometryVariables &Geometry, Node &N, Customise &Custom);
	void readcustom(Customise &Custom);
};

void Input::readgeometry(GeometryVariables &Geometry, Surface &S, Node &N, Customise &Custom)
{
	str line;
	int linenumber(0);
	std::ifstream infile;
	S.setnumber(0);
	N.setnumber(0);
	int col(0);

	infile.open(Custom.getmatfile(), std::ifstream::in);

	if (infile.is_open())
	{
		while (getline(infile, line))
		{

			if (line.empty() || infile.bad())
			{
				linenumber = 0;
				N.setnumber(0);
				getline(infile, line);
			}

			if (infile.good() && line.at(0) != '#')
			{
				line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1"); //remove trailing and leading spaces
				line += " ";

				if(linenumber==0)
				{
					Geometry.setintsurface(S.getnumber(), S);
					S.settype(line);
					Geometry.getintsurface(S.getnumber()).settype(line);
				}
				else if(linenumber==1)
				{
					if (is_number(line))
					{
						S.setmatid(stoi(line));
						Geometry.getintsurface(S.getnumber()).setmatid(stoi(line));
					}
					else
					{
						S.setmatname(line);
						Geometry.getintsurface(S.getnumber()).setmatname(line);
					}
				}
				else
				{
					S.setintnode(N.getnumber(), N);
					Geometry.getintsurface(S.getnumber()).setintnode(N.getnumber(), N);
					str delimiter = " ";
					size_t pos = 0;
					str token;
					col = 0;
					while ((pos = line.find(delimiter)) != str::npos)
					{
						token = line.substr(0, pos);
						assert(col <= 1 && "Error (in readgeometry): Too many columns in input file.\n");

						switch (col)
						{
						case(0):
						{
							assert((std::stoi(token) == 0 || std::stoi(token) == 1) && "Error (in displacement constraints): Constraint must be 0 or 1.\n");
							N.setdx(std::stod(token));
							S.getintnode(N.getnumber()).setdx(std::stod(token));
							Geometry.getintsurface(S.getnumber()).getintnode(N.getnumber()).setdx(std::stod(token));
							break;
						}
						case(1):
						{
							assert((std::stoi(token) == 0 || std::stoi(token) == 1) && "Error (in displacement constraints): Constraint must be 0 or 1.\n");
							N.setdy(std::stod(token));
							S.getintnode(N.getnumber()).setdy(std::stod(token));
							Geometry.getintsurface(S.getnumber()).getintnode(N.getnumber()).setdy(std::stod(token));
							break;
						}
						}
						line.erase(0, pos + delimiter.length());
						col++;
					}
				}
				linenumber++;
				S.setnumber(S.getnumber() + 1);
			}
		}

		infile.close();
	}
	else { std::cerr << "Unable to open geometry file" << std::endl; }

	Geometry.setnumsurfaces(S.getnumber());
	assert(Geometry.getnumsurfaces() > 0 && "Error (in readgeometry): no surfaces assigned.");
};

void Input::readdisplconstraints(GeometryVariables &Geometry, Node &N, Customise &Custom)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int col(0);
	int numdisplconstraints(0);

	infile.open(Custom.getndcfile());
	if (infile.is_open())
	{
		while (getline(infile, line)) 
		{
			if (line.empty())
			{
				getline(infile, line);
			}

			if (infile.good() && line.at(0) != '#')
			{
				line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1"); //remove trailing and leading spaces
				line += " ";

				lines.push_back(line);
				numdisplconstraints++;

				str delimiter = " ";
				size_t pos = 0;
				str token;
				col = 0;
				while ((pos = line.find(delimiter)) != str::npos)
				{

					token = line.substr(0, pos);

					assert(col <= 2 && "Error (in displacement constraints): Too many columns in input file.\n");

					switch (col)
					{
					case(0): 
					{
						N.setnumber(std::stoi(token)); 
						assert(0 <= std::stoi(token) < Geometry.getmaxnumnodes() && "Error (in displacement constraints): Node number invalid.\n");
						break; 
					}
					case(1):
					{
						assert((std::stoi(token) == 0 || std::stoi(token) == 1) && "Error (in displacement constraints): Constraint must be 0 or 1.\n");
						N.setdcx(std::stod(token));
						Geometry.getintnode(N.getnumber()).setdcx(std::stod(token));
						break;
					}
					case(2):
					{
						assert((std::stoi(token) == 0 || std::stoi(token) == 1) && "Error (in displacement constraints): Constraint must be 0 or 1.\n");
						N.setdcy(std::stod(token));
						Geometry.getintnode(N.getnumber()).setdcy(std::stod(token));
						break;
					}
					}
					line.erase(0, pos + delimiter.length());
					col++;
				}
			}

			if (N.getdcx() == 1 || N.getdcy() == 1)
			{
				N.setcflag(1);
				Geometry.getintnode(N.getnumber()).setcflag(1);
			}
		}
		Geometry.setnumconstraints(Geometry.getnumconstraints() + numdisplconstraints);
		assert(Geometry.getnumconstraints() < 0 && "Error (in displacement constraints): Invalid number of constraints.\n");
		infile.close();
	}
	else { std::cerr << "Error (in readdisplacementconstraints): Unable to open file." << std::endl; }
}

void Input::readvelconstraints(GeometryVariables &Geometry, Node &N, Customise &Custom)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int col(0);
	int numvelconstraints(0);

	infile.open(Custom.getnvcfile());
	if (infile.is_open())
	{
		while (getline(infile, line)) 
		{
			if (line.empty())
			{
				getline(infile, line);
			}

			if (infile.good() && line.at(0) != '#')
			{
				line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1"); //remove trailing and leading spaces
				line += " ";

				lines.push_back(line);
				numvelconstraints++;

				str delimiter = " ";
				size_t pos = 0;
				str token;
				col = 0;


				while ((pos = line.find(delimiter)) != str::npos) 
				{

					assert(col <= 2 && "Error (in velocity constraints): Too many columns in input file.\n");

					token = line.substr(0, pos);
					switch (col)
					{
					case(0): 
					{
						assert(0 <= std::stoi(token) < Geometry.getmaxnumnodes() && "Error (in velocity constraints): Node number invalid.\n");
						N.setnumber(std::stoi(token)); 
						break; 
					}
					case(1):
					{
						assert((std::stoi(token) == 0 || std::stoi(token) == 1) && "Error (in velocity constraints): Constraint must be 0 or 1.\n");
						N.setvcx(std::stod(token));
						Geometry.getintnode(N.getnumber()).setvcx(std::stod(token)); 
						break; 
					}
					case(2): 
					{
						assert((std::stoi(token) == 0 || std::stoi(token) == 1) && "Error (in velocity constraints): Constraint must be 0 or 1.\n");
						N.setvcy(std::stod(token));
						Geometry.getintnode(N.getnumber()).setvcy(std::stod(token)); 
						break; 
					}
					}
					line.erase(0, pos + delimiter.length());
					col++;
				}
				if (N.getvcx() == 1 || N.getvcy() == 1)
				{
					N.setcflag(1);
					Geometry.getintnode(N.getnumber()).setcflag(1);
				}
			}
		}
		Geometry.setnumconstraints(Geometry.getnumconstraints() + numvelconstraints);
		assert(Geometry.getnumconstraints() < 0 && "Error (in velocity constraints): Invalid number of constraints.\n");
		infile.close();
	}
	else { std::cerr << "Error (in readvelocityconstraints) Unable to open file." << std::endl; }
}

void Input::readforceconstraints(GeometryVariables &Geometry, Node &N, Customise &Custom)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int col(0);
	int numforceconstraints(0);

	infile.open(Custom.getnfcfile());
	if (infile.is_open())
	{
		while (getline(infile, line)) 
		{
			if (line.empty())
			{
				getline(infile, line);
			}

			if (infile.good() && line.at(0) != '#')
			{
				line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1"); //remove trailing and leading spaces
				line += " ";

				lines.push_back(line);
				numforceconstraints++;

				str delimiter = " ";
				size_t pos = 0;
				str token;
				col = 0;
				while ((pos = line.find(delimiter)) != str::npos)
				{

					assert(col <= 2 && "Error (in force constraints): Too many columns in input file.\n");

					token = line.substr(0, pos);
					switch (col)
					{
					case(0): 
					{
						assert(0 <= std::stoi(token) < Geometry.getmaxnumnodes() && "Error (in force constraints): Node number invalid.\n");
						N.setnumber(std::stoi(token)); 
						Geometry.setintnode(N.getnumber(), N);
						break; 
					}
					case(1): 
					{ 
						assert((std::stoi(token) == 0 || std::stoi(token) == 1) && "Error (in force constraints): Constraint must be 0 or 1.\n");
						N.setfcx(std::stod(token));
						Geometry.getintnode(N.getnumber()).setfcx(std::stod(token)); 
						break; 
					}
					case(2): 
					{ 
						assert((std::stoi(token) == 0 || std::stoi(token) == 1) && "Error (in force constraints): Constraint must be 0 or 1.\n");
						N.setfcy(std::stod(token));
						Geometry.getintnode(N.getnumber()).setfcy(std::stod(token)); 
						break; 
					}
					}
					line.erase(0, pos + delimiter.length());
					col++;
				}
				if (N.getfcx() == 1 || N.getfcy() == 1)
				{
					N.setcflag(1);
					Geometry.getintnode(N.getnumber()).setcflag(1);
				}
			}
		}
		Geometry.setnumconstraints(Geometry.getnumconstraints() + numforceconstraints);
		assert(Geometry.getnumconstraints() < 0 && "Error (in force constraints): Invalid number of constraints.\n");
		infile.close();
	}
	else { std::cerr << "Error (in readforceconstraints): Unable to open file." << std::endl; }
}

void Input::readelements(GeometryVariables &Geometry, Element &E, Customise &Custom)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int col(0);
	int elecounter(0);

	infile.open(Custom.getelementfile());
	if (infile.is_open())
	{
		while (getline(infile, line))
		{
			if (line.empty())
			{
				getline(infile, line);
			}

			if (infile.good() && line.at(0) != '#')
			{
				elecounter++;
			}
		}
		infile.close();
	}
	else { std::cerr << "Error (in readnodes): Unable to open file." << std::endl; }

	Geometry.setnumele(elecounter + 1);
	elecounter = 0;

	infile.open(Custom.getelementfile());
	if (infile.is_open())
	{
		while (getline(infile, line))
		{
			if (line.empty())
			{
				getline(infile, line);
			}

			if (infile.good() && line.at(0) != '#')
			{
				line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1"); //remove trailing and leading spaces
				line += " ";

				for (int i = line.size() - 1; i >= 0; i--)
				{
					if (line[i] == ' '&& line[i] == line[i - 1])
					{
						line.erase(line.begin() + i);
					}
				}

				lines.push_back(line);
				E.setnumber(elecounter);
				Geometry.setintele(E.getnumber(), E);

				str delimiter = " ";
				size_t pos = 0;
				str token;
				col = 0;
				while ((pos = line.find(delimiter)) != str::npos)
				{

					assert(col < Geometry.getmaxnodesperfe() && "Error (in readelements): Max nodes per FE exceeded.\n");

					token = line.substr(0, pos);
					E.setnodenum(col, std::stoi(token));
					Geometry.getintele(E.getnumber()).setnodenum(col, std::stoi(token));
					line.erase(0, pos + delimiter.length());
					col++;
				}
				E.setnumnodes(col);
				Geometry.getintele(E.getnumber()).setnumnodes(E.getnumnodes());

				elecounter++;
			}
		}
		infile.close();
	}
	else { std::cerr << "Error (in readelements): Unable to open file." << std::endl; }

	Geometry.setelementdata(lines);
	//Geometry.setnumele(lines.size());
	Geometry.setmaxnumele(10 * Geometry.getnumele());
	assert(Geometry.getmaxnumele() >= 0 && "Error (in readelements): Maximum number of elements invalid.\n");
	if (Geometry.getmaxnumele() < 10)
	{
		std::cerr << "Warning (in readelements): Maximum number of elements may be exceeded.\n";
	}
}

void Input::readnodes(GeometryVariables &Geometry, Node &N, Customise &Custom, PropertyVariables &Property)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int col(0);
	N.setnumber(0);
	Geometry.setnumboundary(0);
	int nodecounter(0);

	Geometry.setnodesperfe(3);
	Geometry.setmaxnodesperfe(4);

	infile.open(Custom.getnodefile());
	if (infile.is_open())
	{
		while (getline(infile, line))
		{
			if (line.empty())
			{
				getline(infile, line);
			}

			if (infile.good() && line.at(0) != '#')
			{
				nodecounter++;
			}
		}
		infile.close();
	}
	else { std::cerr << "Error (in readnodes): Unable to open file." << std::endl; }

	assert(nodecounter >= 0 && "Error (in readnodes): Number of nodes invalid.");
	if (nodecounter < 10)
	{
		std::cerr << "Warning (in readnodes): Very few nodes set.\n";
	}

	Geometry.setnumnodes(nodecounter);
	nodecounter = 0;

	infile.open(Custom.getnodefile());
	if (infile.is_open())
	{
		while (getline(infile, line))
		{
			if (line.empty())
			{
				getline(infile, line);
			}

			if (infile.good() && line.at(0) != '#')
			{
				line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1"); //remove trailing and leading spaces
				line += " ";

				for (int i = line.size() - 1; i >= 0; i--)
				{
					if (line[i] == ' '&& line[i] == line[i - 1])
					{
						line.erase(line.begin() + i);
					}
				}

				lines.push_back(line);
				N.setnumber(nodecounter);
				Geometry.setnumboundary(Geometry.getnumboundary() + 1);

				//set boundary flag to 1 by default
				N.setboundaryflag(1);
				Geometry.getintnode(N.getnumber()).setboundaryflag(1);

				str delimiter = " ";
				size_t pos = 0;
				double token;
				col = 0;

				while ((pos = line.find(delimiter)) != str::npos)
				{
					token = std::stod(line.substr(0, pos));

					assert(col > 3 && "Error (in readnodes): Too many columns in input file.\n");

					switch (col)
					{
					case(0):
					{
						N.setdx(token);
						Geometry.getintnode(N.getnumber()).setdx(token);
						break;
					}
					case(1):
					{
						N.setdy(token); Geometry.getintnode(N.getnumber()).setdy(token);
						break;
					}
					case(2):
					{
						assert(token < Property.getnummats() && "Error (in readnodes): material id invalid.\n");
						N.setmatid(int(token)); 
						Geometry.getintnode(N.getnumber()).setmatid(int(token));
						break;
					}
					case(3):
					{
						assert((int(token) == 0 || int(token) == 1) && "Error (in readnodes): boundary flag must be either 0 or 1.\n");
						N.setboundaryflag(int(token)); 
						Geometry.getintnode(N.getnumber()).setboundaryflag(int(token));
						break;
					}

					line.erase(0, pos + delimiter.length());
					col++;
					}

					Geometry.setintnode(N.getnumber(), N);
					nodecounter++;
				}
			}
		}
		Geometry.setnodedata(lines);
		Geometry.setmaxnumnodes(50 * Geometry.getnumnodes());
		assert(Geometry.getmaxnumnodes() > 0 && "Error (in readnodes): Maximum number of nodes invalid.\n");
		infile.close();
	}
	else
	{
		std::cerr << "Error (in readnodes): Unable to open file" << std::endl;
	}
}

void Input::readnodalforces(GeometryVariables &Geometry, Node &N, Customise &Custom)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int col(0);
	int numinitforces(0);

	infile.open(Custom.getnodeffile());
	if (infile.is_open())
	{
		while (getline(infile, line))
		{
			if (line.empty())
			{
				getline(infile, line);
			}

			if (infile.good() && line.at(0) != '#')
			{
				line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1"); //remove trailing and leading spaces
				line += " ";

				for (int i = line.size() - 1; i >= 0; i--)
				{
					if (line[i] == ' '&& line[i] == line[i - 1])
					{
						line.erase(line.begin() + i);
					}
				}

				lines.push_back(line);
				numinitforces++;

				str delimiter = " ";
				size_t pos = 0;
				str token;
				col = 0;
				while ((pos = line.find(delimiter)) != str::npos)
				{
					if (col > 2)
					{
						std::cerr << "Error (in nodal forces): Too many columns in input file." << std::endl;
					}

					token = line.substr(0, pos);
					switch (col)
					{
					case(0): { N.setnumber(std::stoi(token)); break; }
					case(1): { N.setfx(std::stod(token)); Geometry.getintnode(N.getnumber()).setfx(std::stod(token)); break; }
					case(2): { N.setfy(std::stod(token)); Geometry.getintnode(N.getnumber()).setfy(std::stod(token)); break; }
					}
					line.erase(0, pos + delimiter.length());
					col++;
				}
				if (N.getfx() != 0 || N.getfy() != 0)
				{
					N.setinitflag(1);
					Geometry.getintnode(N.getnumber()).setinitflag(1);
				}
			}
		}
		Geometry.setnuminit(Geometry.getnuminit() + numinitforces);
		infile.close();
	}
	else { std::cerr << "Error (in readforces): Unable to open file." << std::endl; }
}

void Input::readnodalvelocities(GeometryVariables &Geometry, Node &N, Customise &Custom)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int col(0);
	int numinitvel(0);

	infile.open(Custom.getnodevfile());
	if (infile.is_open())
	{
		while (getline(infile, line))
		{
			if (line.empty())
			{
				getline(infile, line);
			}

			if (infile.good() && line.at(0) != '#')
			{
				line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1"); //remove trailing and leading spaces
				line += " ";

				for (int i = line.size() - 1; i >= 0; i--)
				{
					if (line[i] == ' '&& line[i] == line[i - 1])
					{
						line.erase(line.begin() + i);
					}
				}

				lines.push_back(line);
				numinitvel++;

				str delimiter = " ";
				size_t pos = 0;
				str token;
				col = 0;

				while ((pos = line.find(delimiter)) != str::npos)
				{
					if (col > 2)
					{
						std::cerr << "Error (in nodal velocities): Too many columns in input file." << std::endl;
					}

					token = line.substr(0, pos);
					switch (col)
					{
					case(0): {N.setnumber(std::stoi(token)); break; }
					case(1): { Geometry.getintnode(N.getnumber()).setvx(std::stod(token)); break; }
					case(2): {Geometry.getintnode(N.getnumber()).setvy(std::stod(token)); break; }
					}
					line.erase(0, pos + delimiter.length());
					col++;
				}
				if (N.getvx() != 0 || N.getvy() != 0)
				{
					N.setinitflag(1);
					Geometry.getintnode(N.getnumber()).setinitflag(1);
				}
			}
		}
		Geometry.setnuminit(Geometry.getnuminit() + numinitvel);
		infile.close();
	}
	else { std::cerr << "Error (in nodal velocities): Unable to open file." << std::endl; }
}

void Input::readuser(UserVariables &User, Customise &Custom)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int linenumber(0);

	infile.open(Custom.getuserfile(), std::ifstream::in);
	if (infile.is_open()) 
	{
		while (getline(infile, line)) 
		{
			if (line.empty())
			{
				getline(infile, line);
			}

			if (infile.good() && line.at(0) != '#')
			{
				line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1"); //remove trailing and leading spaces
				line += " ";

				for (int i = line.size() - 1; i >= 0; i--)
				{
					if (line[i] == ' '&& line[i] == line[i - 1])
					{
						line.erase(line.begin() + i);
					}
				}

				lines.push_back(line);
				switch (linenumber)
				{
				case(0):{User.setmakenewfile(line); break; }
				case(1): {User.setoutfile(line); break; }
				case(2):{User.setlogin(line); break; }
				case(3):{User.sethpcdir(line); break; }
				case(4):{User.setlocaldir(line); break; }
				case(5):{User.setqsubdir(line); break; }
				case(6): {User.setwindowsmayavidir(line); break; }
				case(7): {User.setbashmayavidir(line); break; }
				case(8): {User.setrunlocal(line); break; }
				case(9): {User.setwindowsgiddir(line); break; }
				case(10): {User.setbashgiddir(line); break; }
				case(11): {User.setwindowsydir(line); break; }
				case(12): {User.setbashydir(line); break; }
				case(13): {User.setwindowsparaviewdir(line); break; }
				case(14): {User.setbashparaviewdir(line); break; }
				}
				linenumber++;
			}
		}

		User.setuserdata(lines);
		infile.close();
	}
	else { std::cerr << "Error (in readuser): Unable to open file." << std::endl; }
}

void Input::readcustom(Customise &Custom)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int linenumber(0);

	infile.open("../user/custom.txt", std::ifstream::in);
	if (infile.is_open())
	{
		while (getline(infile, line))
		{
			if (line.empty())
			{
				getline(infile, line);
			}

			if (infile.good() && line.at(0) != '#')
			{
				line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1"); //remove trailing and leading spaces
				line += " ";

				for (int i = line.size() - 1; i >= 0; i--)
				{
					if (line[i] == ' '&& line[i] == line[i - 1]) //added equal sign
					{
						line.erase(line.begin() + i);
					}
				}

				switch (linenumber)
				{
				case(0): {Custom.setuserfile(line); break; }
				case(1): {Custom.setmatfile(line); break; }
				case(2): {Custom.setctrlfile(line); break; }	
				case(3): {Custom.setnodefile(line); break; }	
				case(4): {Custom.setnodevfile(line); break; }	
				case(5): {Custom.setnodeffile(line); break; }	
				case(6): {Custom.setndcfile(line); break; }	
				case(7): {Custom.setnvcfile(line); break; }	
				case(8): {Custom.setnfcfile(line); break; }	
				case(9): {Custom.setinitfile(line); break; }	
				case(10): {Custom.setelementfile(line); break; }
				case(11): {Custom.setgeometryfile(line); break; }
				}
				lines.push_back(line);
				linenumber++;
			}
		}

		Custom.setdata(lines);
		infile.close();
	}
	else { std::cerr << "Error (in readcustom): Unable to open custom.txt ." << std::endl; }
}

void Input::readproblem(ProblemVariables &Problem, PropertyVariables &Property, Customise &Custom, GeometryVariables &Geometry)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int linenumber(0);

	infile.open(Custom.getctrlfile(), std::ifstream::in);
	if (infile.is_open())
	{
		while (getline(infile, line)) 
		{
			if (line.empty())
			{
				getline(infile, line);
			}

			if (infile.good() && line.at(0) != '#')
			{
				line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1"); //remove trailing and leading spaces
				line += " ";

				switch (linenumber)
				{
				case(0):
				{
					Problem.setrealsimtime(std::stod(line));
					assert(Problem.getrealsimtime() > 0 && "Error (in readproblem): Real simulation time invalid.\n");
					break;
				}
				case(1): 
				{
					Problem.setgravity(std::stod(line)); 
					if (int(Problem.getgravity()) == 0) { std::cerr << "Warning (in readproblem): Small gravity.\n" << std::endl; }
					break; 
				}
				case(2): 
				{
					Problem.setnumberoutfiles(stoi(line));
					assert(Problem.getnumberoutfiles() >= 0 && "Error (in readproblem): Number of outfiles invalid.\n");
					Problem.setoutfreq(int((double) Problem.getmaxndt() / (double) Problem.getnumberoutfiles())); 
					assert(Problem.getoutfreq() >= 0 && "Error (in readproblem): Output frequency invalid.\n");
					break; 
				}
				case(3): 
				{
					Problem.setrestartfreq(Problem.getmaxndt());
					break; 
				}
				case(4): 
				{
					Problem.setgravitysettlingstage(std::stoi(line)); 
					break; 
				}
				case(5): 
				{
					Problem.setloadrampingstage(std::stoi(line)); 
					break; 
				}
				case(6): 
				{
					Problem.setmaxdim(std::stoi(line)); 
					assert(4 > Problem.getmaxdim() > 0 && "Error (in readproblem): Invalid maximum dimension.\n");
					break; 
				}
				case(7): 
				{
					Problem.setmaxforce(std::stod(line));
					if (Problem.getmaxforce() < 100) { std::cerr << "Warning (in readproblem): Maximum force may be exceeded.\n"; }
					break; 
				}
				case(8): 
				{
					Problem.setmaxv(std::stoi(line)); 
					if (Problem.getmaxv() < 1) { std::cerr << "Warning (in readproblem): Maximum velocity may be exceeded.\n"; }
					break; 
				}
				case(9): 
				{
					Problem.setmaxstress(std::stod(line)); 
					if (Problem.getmaxstress() < 100) { std::cerr << "Warning (in readproblem): Maximum stress may be exceeded.\n"; }
					break; 
				}
				case(10): 
				{
					Problem.setmaxdispl(std::stoi(line)); 
					if (Problem.getmaxstress() < 1) { std::cerr << "Warning (in readproblem): Maximum displacement may be exceeded.\n"; }
					break; 
				}
				case(11): 
				{
					Problem.setmaxjointaperture(std::stod(line)); 
					break; 
				}
				case(12): 
				{
					Problem.setmaxcontcouples(std::stoi(line)); 
					if (Problem.getmaxcontcouples() < 10) { std::cerr << "Warning (in readproblem): Maximum number of joint couples may be exceeded.\n"; }
					break; 
				}
				case(13): {Problem.setaccuracy(line); break; }
				case(14): {Problem.setfricmodel(line); break; }
				case(15): {Problem.setinitialaperturecorr(line); break; }
				case(16): {Problem.setgenoutonly(line); break; }
				}

				lines.push_back(line);
				linenumber++;
			}
		}
		infile.close();
	}
	else { std::cerr << "Error (in readproblem): Unable to open file." << std::endl; }

	assert(Problem.getcurrentdt() > -1 && "Error (in readproblem): Current time steps too invalid.\n");

	Problem.setproblemdata(lines); // entire data

	Problem.setdt(std::floor(Property.getdtcrit()));
	Problem.setmaxndt(std::ceil(Problem.getrealsimtime() / Problem.getdt())); //number of time steps

	std::cout << "Time step: " << Problem.getdt() << std::endl;
	std::cout << "Real simulation time: " << Problem.getrealsimtime() << std::endl;
	std::cout << "Number of time steps: " << Problem.getrealsimtime() / Problem.getdt()
		<< ", adjusted to: " << Problem.getmaxndt() << std::endl;

	assert(int(Problem.getmaxndt()) >= 0 && "Error (in readproblem): Maximum number of time steps invalid.\n");

	Problem.setcurrentndt(0); //current time step
	Problem.setcoordsize(1);
	Problem.setcurrentt((double)0); //current time
	Problem.setcontcouples(0); //current number of contacting couples

	Problem.setbuffersize(0.1 * Geometry.getminedge()); //buffer size
	assert(Problem.getbuffersize() > 0 && "Error (in readproblem): Invalid buffer size.\n");
}

void Input::readmaterial(Material Mat, PropertyVariables &Property, ProblemVariables &Problem, Customise &Custom, GeometryVariables &Geometry)
{
	str line;
	int linenumber(0);
	std::ifstream infile;
	Mat.setnumber(0);
	bool resetmat = false;

	infile.open(Custom.getmatfile(), std::ifstream::in);

	if (infile.is_open())
	{
		while (getline(infile, line))
		{

			if (line.empty() || infile.bad())
			{
				linenumber = 0;
				getline(infile, line);
			}

			if (infile.good() && line.at(0) != '#')
			{
				line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1"); //remove trailing and leading spaces
				line += " ";


				switch (linenumber)
				{
				case(0): 
				{ 
					Mat.setname(line); 
					Property.setmatname(Mat.getnumber(), Mat.getname());
					Property.setintmat(Mat.getnumber(), Mat);
					Property.setstrmat(Mat.getname(), Mat);

					Mat.setviscousdamping(0); //????
					Mat.setloadamplitude(0);
					Mat.setamplitude(0);
					Mat.setamplitudeinc(0);

					Property.getintmat(Mat.getnumber()).setviscousdamping(0);
					Property.getstrmat(Mat.getname()).setviscousdamping(0);

					Property.getintmat(Mat.getnumber()).setloadamplitude(0);
					Property.getstrmat(Mat.getname()).setloadamplitude(0);

					Property.getintmat(Mat.getnumber()).setamplitude(0);
					Property.getstrmat(Mat.getname()).setamplitude(0);

					Property.getintmat(Mat.getnumber()).setamplitudeinc(0);
					Property.getstrmat(Mat.getname()).setamplitudeinc(0);
					break;
				}

				case(1): 
				{ 
					Mat.setdensity(std::stod(line)); 
					Property.getintmat(Mat.getnumber()).setdensity(std::stod(line));
					Property.getstrmat(Mat.getname()).setdensity(std::stod(line));
					assert(Mat.getdensity() > 0 && "Error (in readmaterial): Invalid density.\n");
					break; 
				}
				case(2): 
				{ 
					Mat.setyoungsmod(std::stod(line)); 
					Property.getintmat(Mat.getnumber()).setyoungsmod(std::stod(line));
					Property.getstrmat(Mat.getname()).setyoungsmod(std::stod(line));
					assert(Mat.getdensity() > 0 && "Error (in readmaterial): Invalid young's modulus.\n");

					Mat.setelasticpenalty(10*Mat.getyoungsmod());
					Property.getintmat(Mat.getnumber()).setelasticpenalty(10 * Mat.getyoungsmod());
					Property.getstrmat(Mat.getname()).setelasticpenalty(10 * Mat.getyoungsmod());
					assert(Mat.getelasticpenalty() > 0 && "Error (in readmaterial): Invalid elastic penalty.\n");

					Mat.setcontactpenalty(10 * Mat.getyoungsmod());
					Property.getintmat(Mat.getnumber()).setcontactpenalty(10 * Mat.getyoungsmod());
					Property.getstrmat(Mat.getname()).setcontactpenalty(10 * Mat.getyoungsmod());
					assert(Mat.getcontactpenalty() > 0 && "Error (in readmaterial): Invalid contact penalty.\n");

					Mat.setdtcrit(Geometry.getminelev() * Mat.getdensity() / Mat.getelasticpenalty());
					Property.getintmat(Mat.getnumber()).setdtcrit(Geometry.getminelev() * Mat.getdensity() / Mat.getelasticpenalty());
					Property.getstrmat(Mat.getname()).setdtcrit(Geometry.getminelev() * Mat.getdensity() / Mat.getelasticpenalty());
					assert(Mat.getdtcrit() > 0 && "Error (in readmaterial): Invalid critical timestep.\n");
					break; 
				}
				case(3): 
				{ 
					Mat.setnu(std::stod(line)); 
					Property.getintmat(Mat.getnumber()).setnu(std::stod(line));
					Property.getstrmat(Mat.getname()).setnu(std::stod(line));
					assert(Mat.getnu() > 0 && "Error (in readmaterial): Invalid Poisson's ratio.\n");

					Mat.setlambda(Mat.getnu() * Mat.getyoungsmod() / (1 + Mat.getnu()) / (1 - 2 * Mat.getnu()));
					Property.getintmat(Mat.getnumber()).setlambda(Mat.getnu() * Mat.getyoungsmod() / (1 + Mat.getnu()) / (1 - 2 * Mat.getnu()));
					Property.getstrmat(Mat.getname()).setlambda(Mat.getnu() * Mat.getyoungsmod() / (1 + Mat.getnu()) / (1 - 2 * Mat.getnu()));
					assert(Mat.getlambda() > 0 && "Error (in readmaterial): Invalid lame constant lambda.\n");

					Mat.setmu(Mat.getyoungsmod() / (2 * (1 + Mat.getnu())));
					Property.getintmat(Mat.getnumber()).setmu(Mat.getyoungsmod() / (2 * (1 + Mat.getnu())));
					Property.getstrmat(Mat.getname()).setmu(Mat.getyoungsmod() / (2 * (1 + Mat.getnu())));
					assert(Mat.getmu() > 0 && "Error (in readmaterial): Invalid lame constant mu.\n");

					Mat.setmassdamping(sqrt(Mat.getyoungsmod()*Mat.getdensity()) * Geometry.getminedge());
					Property.getintmat(Mat.getnumber()).setmassdamping(sqrt(Mat.getyoungsmod()*Mat.getdensity()) * Geometry.getminedge());
					Property.getstrmat(Mat.getname()).setmassdamping(sqrt(Mat.getyoungsmod()*Mat.getdensity()) * Geometry.getminedge());
					assert(Mat.getmassdamping() > 0 && "Error (in readmaterial): Invalid mass damping.\n");
					break; 
				}
				case(6): 
				{ 
					Mat.setmodeienergyrate(std::stod(line)); 
					Property.getintmat(Mat.getnumber()).setmodeienergyrate(std::stod(line));
					Property.getstrmat(Mat.getname()).setmodeienergyrate(std::stod(line));
					assert(Mat.getmodeienergyrate() > 0 && "Error (in readmaterial): Invalid mode 1 energy rate.\n");
					break; 
				}
				case(7): 
				{ 
					Mat.setmodeiienergyrate(std::stod(line)); 
					Property.getintmat(Mat.getnumber()).setmodeiienergyrate(std::stod(line));
					Property.getstrmat(Mat.getname()).setmodeiienergyrate(std::stod(line));
					assert(Mat.getmodeiienergyrate() > 0 && "Error (in readmaterial): Invalid mode 2 energy rate.\n");
					break; 
				}
				case(8):
				{ 
					Mat.settensilestrength(std::stod(line)); 
					Property.getintmat(Mat.getnumber()).settensilestrength(std::stod(line));
					Property.getstrmat(Mat.getname()).settensilestrength(std::stod(line));
					assert(Mat.gettensilestrength() > 0 && "Error (in readmaterial): Invalid tensile strength.\n");
					break; 
				}
				case(9): 
				{ 
					Mat.setinternalfriction(std::stod(line)); 
					Property.getintmat(Mat.getnumber()).setinternalfriction(std::stod(line));
					Property.getstrmat(Mat.getname()).setinternalfriction(std::stod(line));
					assert(Mat.getinternalfriction() > 0 && "Error (in readmaterial): Invalid internal friction.\n");
					break; 
				}
				case(10): 
				{ 
					Mat.setinternalcohesion(std::stod(line)); 
					Property.getintmat(Mat.getnumber()).setinternalcohesion(std::stod(line));
					Property.getstrmat(Mat.getname()).setinternalcohesion(std::stod(line));
					assert(Mat.getinternalcohesion() > 0 && "Error (in readmaterial): Invalid internal cohesion.\n");
					break; 
				}
				case(11): 
				{ 
					Mat.setporefluidpressure(std::stod(line));
					Property.getintmat(Mat.getnumber()).setporefluidpressure(std::stod(line));
					Property.getstrmat(Mat.getname()).setporefluidpressure(std::stod(line));
					break; 
				}
				case(12): 
				{ 
					Mat.setjointfriction(std::stod(line)); 
					Property.getintmat(Mat.getnumber()).setjointfriction(std::stod(line));
					Property.getstrmat(Mat.getname()).setjointfriction(std::stod(line));
					break; 
				}
				case(13): 
				{ 
					Mat.setjrc0(std::stod(line)); 
					Property.getintmat(Mat.getnumber()).setjrc0(std::stod(line));
					Property.getstrmat(Mat.getname()).setjrc0(std::stod(line));
					break; 
				}
				case(14): 
				{ 
					Mat.setjcs0(std::stod(line)); 
					Property.getintmat(Mat.getnumber()).setjcs0(std::stod(line));
					Property.getstrmat(Mat.getname()).setjcs0(std::stod(line));
					break; 
				}
				case(15): 
				{ 
					Mat.setjointsamplesize(std::stod(line)); 
					Property.getintmat(Mat.getnumber()).setjointsamplesize(std::stod(line));
					Property.getstrmat(Mat.getname()).setjointsamplesize(std::stod(line));
					break; 
				}
				case(16): 
				{ 
					Mat.setinterfacefriction(std::stod(line)); 
					Property.getintmat(Mat.getnumber()).setinterfacefriction(std::stod(line));
					Property.getstrmat(Mat.getname()).setinterfacefriction(std::stod(line));
					assert(Mat.getinterfacefriction() > 0 && "Error (in readmaterial): Invalid interface friction.\n");
					break; 
				}
				case(17): 
				{ 
					Mat.setproblems2d(line); 
					Property.getintmat(Mat.getnumber()).setproblems2d(line);
					Property.getstrmat(Mat.getname()).setproblems2d(line);
					break; 
				}
				case(18): 
				{ 
					Mat.setfractureflag(stoi(line)); 
					Property.getintmat(Mat.getnumber()).setfractureflag(stoi(line));
					Property.getstrmat(Mat.getname()).setfractureflag(stoi(line));
					break; 
				}
				case(19): 
				{ 
					Mat.setnumberofmeshrefinements(stoi(line)); 
					Property.getintmat(Mat.getnumber()).setnumberofmeshrefinements(stoi(line));
					Property.getstrmat(Mat.getname()).setnumberofmeshrefinements(stoi(line));
					break; 
				}
				case(20): 
				{ 
					Mat.setpropertytype(stoi(line)); 
					Property.getintmat(Mat.getnumber()).setpropertytype(stoi(line));
					Property.getstrmat(Mat.getname()).setpropertytype(stoi(line));
					break; 
				}
				case(21): 
				{ 
					Mat.setslidingfriction(stod(line)); 
					Property.getintmat(Mat.getnumber()).setslidingfriction(stod(line));
					Property.getstrmat(Mat.getname()).setslidingfriction(stod(line));
					break; 
				}
				case(22): 
				{ 
					Mat.setjointproperty(stoi(line)); 
					Property.getintmat(Mat.getnumber()).setjointproperty(stoi(line));
					Property.getstrmat(Mat.getname()).setjointproperty(stoi(line));

					//store the material in the property list
					break; 
				}
				}
				linenumber++;
				Mat.setnumber(Mat.getnumber() + 1);
			}
		}

		infile.close();
	}
	else { std::cerr << "Unable to open material file" << std::endl; }

	assert(Mat.getnumber() > 0 && "Error (in readmaterial): No materials assigned.\n");

	Property.setnummats(Mat.getnumber());
	Property.setnumprop(Property.getnummats());
	Property.setmaxnumprop(Property.getnummats()+1);
};

class Output
{
	friend class Material;
	friend class PropertyVariables;

private:
	std::ofstream file;
	str yname; //set name of the output .y file
	str geoname; //set name of the output .geo file
	str meshname; //set name for gmsh mesh

public:

	Output() {};
	~Output() {};

	str getyname() { return this->yname; }
	void setyname(str yname) { this->yname = yname; }

	str getgeoname() { return this->geoname; }
	void setgeoname(str geoname) { this->geoname = geoname; }

	str getmeshname() { return this->meshname; }
	void setmeshname(str meshname) { this->meshname = meshname; }

	void writegeo(int argc, char **argv);
	void writecontrol(ProblemVariables &Problem);
	void writeelements(GeometryVariables &Geometry);
	void writeinteractions(ProblemVariables &Problem);
	void writenodes(GeometryVariables &Geometry);
	void writeproperties(PropertyVariables &Property);
	void writeoutput(PropertyVariables &Property, GeometryVariables &Geometry, Element &Ele, Node &N, ProblemVariables &Problem);
};		

void Output::writegeo(int argc, char **argv)
{
	gmsh::initialize(argc, argv);

	gmsh::model::add(this->getgeoname());

	//Points
	gmsh::model::geo::addPoint(0, 0, 0, 0.1, 1);
	gmsh::model::geo::addPoint(1, 0, 0, 0.1, 2);
	gmsh::model::geo::addPoint(1, 1, 0, 0.1, 3);
	gmsh::model::geo::addPoint(0, 1, 0, 0.1, 4);

	//Lines
	gmsh::model::geo::addLine(1, 2, 1);
	gmsh::model::geo::addLine(2, 3, 2);
	gmsh::model::geo::addLine(3, 4, 3);

	// try automatic assignement of tag

	int line4 = gmsh::model::geo::addLine(4, 1);
	gmsh::model::geo::addCurveLoop({ 1, 2, 3, line4 }, 1);
	gmsh::model::geo::addPlaneSurface({ 1 }, 6);
	gmsh::model::geo::synchronize();
	gmsh::model::mesh::generate(2);
	gmsh::write(this->getmeshname());
}

void Output::writecontrol(ProblemVariables &Problem)
{
	this->file.open(this->getyname(), std::ios::out | std::ios::trunc | std::ios::binary); ;
	if (this->file.is_open())
	{
		this->file << "/* Control */\n";
		this->file << "/YD/YDC/MCSTEP "		 						<< int(Problem.getmaxndt())				<< "\n"; //Maximum number of steps
		this->file << "/YD/YDC/NCSTEP "								<< Problem.getcurrentdt()				<< "\n"; //Current / Actual number of time steps (cannot be greater than /YD/YDC/MCSTEP)
		this->file << "/YD/YDC/DCGRAY "			<< std::scientific  << Problem.getgravity()					<< "\n"; //Acceleration of gravity (in y direction)
		this->file << "/YD/YDC/DCSIZC "			<< std::scientific	<< Problem.getcoordsize()				<< "\n"; //Maximum size of coordinate in any direction (size of physical space – corresponding outputs are normalized using this value)
		this->file << "/YD/YDC/DCSIZF "			<< std::scientific	<< Problem.getmaxforce()				<< "\n"; //Maximum size of force in any direction(corresponding outputs are normalized using this value)
		this->file << "/YD/YDC/DCSIZS "			<< std::scientific	<< Problem.getmaxstress()				<< "\n"; //Maximum stress (buffer size)
		this->file << "/YD/YDC/DCSIZV "			<< std::scientific	<< Problem.getmaxv()					<< "\n"; //Maximum size of velocity in any direction (corresponding outputs are normalized using this value)
		this->file << "/YD/YDC/DCSIZD "			<< std::scientific	<< Problem.getmaxdispl()				<< "\n"; //Maximum size of displacement in any direction
		this->file << "/YD/YDC/DCSIZA "			<< std::scientific	<< Problem.getmaxjointaperture()		<< "\n"; //Maximum joint aperture
		this->file << "/YD/YDC/DCSTEC "			<< std::scientific	<< Problem.getdt()						<< "\n"; //Size of the time step (see in the book how to calculate critical-maximum time step)
		this->file << "/YD/YDC/DCTIME "			<< std::scientific	<< Problem.getcurrentt()				<< "\n"; //Current time. i.e. time at start of this run.
		this->file << "/YD/YDC/DCRMPT "			<< std::scientific	<< Problem.getloadrampingstage()		<< "\n"; //Load ramping stage
		this->file << "/YD/YDC/DCGRST "			<< std::scientific  << Problem.getgravitysettlingstage()	<< "\n"; //Gravity settling stage
		this->file << "/YD/YDC/ICOUTF "								<< int(Problem.getoutfreq())			<< "\n"; //Output frequency – every so many time steps complete state of the system is recorded in a file with extension .ym which can be visualized using M program, which is FEM/DEM Visualizer accompanying Y program.
		this->file << "/YD/YDC/ICOUTI "								<< Problem.getcurrentit()				<< "\n"; //Current number of iterations; Serial number of first output associated with this run.
		this->file << "/YD/YDC/ICSAVF "								<< int(Problem.getrestartfreq())		<< "\n"; //Restart save frequency
		this->file << "/YD/YDC/ICOUTP "								<< Problem.getnumchars()				<< "\n"; //Number of characters for each number in output file (for example: three characters is equivalent to six significant digits)
		this->file << "/YD/YDC/ICFMTY "								<< 0									<< "\n"; //
		this->file << "/YD/YDC/ICIATY "								<< 0									<< "\n"; //
		this->file.close();
	}
	else { std::cerr << "In writecontrol: Unable to open " << this->getyname() << std::endl; }
}

void Output::writeelements(GeometryVariables &Geometry)
{
	this->file.open(this->getyname(), std::ios::out | std::ios::app | std::ios::binary); ;
	if (this->file.is_open())
	{
		this->file << "/* Elements */\n";
		this->file << "/YD/YDE/MELEM "  << Geometry.getmaxnumele()									<< "\n"; //Maximum number of finite elements(with fracture the actual number of finite elements increases during the run, however, it should not exceed this number)
		this->file << "/YD/YDE/NELEM "  << Geometry.getnumele()										<< "\n"; //Actual number of finite elements at the beginning of this run.
		this->file << "/YD/YDE/MELST "  << Geometry.getmaxstatevar()								<< "\n"; //Maximum number of state variables per finite element.
		this->file << "/YD/YDE/NELST "  << Geometry.getstatevar()									<< "\n"; //Actual number of state variables per finite element.
		this->file << "/YD/YDE/MELNO "  << Geometry.getmaxnodesperfe()								<< "\n"; //Maximum number of nodes per finite element.
		this->file << "/YD/YDE/NELNO "  << Geometry.getnodesperfe()									<< "\n"; //Actual number of nodes per finite element.
		this->file << "/YD/YDE/D2ELST " << 21 << " " << Geometry.getnumele() << 0					<< "\n"; //[MELST][MELEM] array containing state variables for all finite elements.
		
		this->file << "/YD/YDE/I1ELCF " << Geometry.getnumele()										<< "\n"; //Head of a list of contacting couples for each finite element.
		for (int i = 0; i != Geometry.getnumele(); i++) { this->file << "-1 "; }
		this->file << "\n";

		this->file << "/YD/YDE/I1ELTY " << Geometry.getnumele()										<< "\n"; //Head of a list of contacting couples for each finite element.
		for (int i = 0; i != Geometry.getnumele(); i++) { this->file << "-1 "; }
		this->file << "\n";

		this->file << "/YD/YDE/I1ELPR " << Geometry.getnumele()										<< "\n"; //Set of properties associated with each element.
		for (int i = 0; i != Geometry.getnumele(); i++) { this->file << "0 "; }
		this->file << "\n";

		this->file << "/YD/YDE/I2ELTO " << 21 << " " << Geometry.getnodesperfe() << " " << Geometry.getnumele(); //[MELNO][MELEM] topology array containing nodes for each finite element.
		for (int i = 0; i != Geometry.getnumele(); i++) //get the ith element
		{ 
			for (int j = 0; j != Geometry.getintele(i).getnumnodes(); j++)
			{
				this->file << Geometry.getintele(i).getnodenum(j) << " ";
			}
			this->file << "\n";
		}

		this->file << "/YD/YDE/I2ELJP " << 21 << " " << Geometry.getnodesperfe() << " " << Geometry.getnumele()	<< "\n";
		this->file.close();
	}
	else { std::cerr << "Error (in writeelements): Unable to open " << this->getyname() << std::endl; }
}

void Output::writeinteractions(ProblemVariables &Problem)
{
	this->file.open(this->getyname(), std::ios::out | std::ios::app | std::ios::binary); ;
	if (this->file.is_open())
	{
		this->file << "/* Interactions */\n";
		this->file << "/YD/YDI/MICOUP "						<< Problem.getmaxcontcouples()	<< "\n"; //Maximum number of contacting couples of finite elements.
		this->file << "/YD/YDI/NICOUP "						<< Problem.getcontcouples()		<< "\n"; //Actual number of contacting couples of finite elements (always set to zero)
		this->file << "/YD/YDI/IIECFF "						<< -2							<< "\n"; //Internal variable used for contact (always set to -2)
		this->file << "/YD/YDI/DIEDI "	<< std::scientific	<< 200							<< "\n"; //Internal variable (travel since last detection) used to trigger contact detection (always set to a large number in order to trigger contact detection immediately at the start of this run)
		this->file << "/YD/YDI/DIEZON " << std::scientific  << Problem.getbuffersize()		<< "\n"; //Contact detection buffer size, Size of the buffer around each finite element for contact detection purposes (usually 1/5 of the size of the smallest finite element)
		this->file << "/YD/YDI/D1IESL "						<< 0							<< "\n"; //[MELEM] array used to store sliding distance between couples in contact
		this->file << "/YD/YDI/I1IECN "						<< 0							<< "\n"; //[MICOUP] array used to store next contacting couple in the list.
		this->file << "/YD/YDI/I1IECT "						<< 0							<< "\n"; //[MICOUP] array used to store target finite element for each contacting couple.
		this->file.close();
	}
	else
	{ std::cerr << "Error (in writeinteractions): Unable to open " << this->getyname() << std::endl; }
}

void Output::writenodes(GeometryVariables &Geometry)
{
	this->file << "/* Nodes */\n";
	this->file.open(this->getyname(), std::ios::out | std::ios::app | std::ios::binary); ;
	if (this->file.is_open())
	{
	this->file << "/YD/YDN/MNODIM "		<< Geometry.getmaxdim()												<< "\n";	//Maximum number of degrees of freedom per node (usually 2 in 2D)
	this->file << "/YD/YDN/NNODIM "		<< Geometry.getdim()												<< "\n";	//Actual number of degrees of freedom per node (usually 2 in 2D)
	this->file << "/YD/YDN/MNOPO "		<< Geometry.getmaxnumnodes()										<< "\n";	//Maximum number of nodes.
	this->file << "/YD/YDN/NNOPO "		<< Geometry.getnumnodes()											<< "\n";	//Actual number of nodes.
	
	this->file << "/YD/YDN/D2NCC "		<< 21 << " "<< Geometry.getdim() << " " << Geometry.getnumnodes()	<< "\n"; //[MNODIM][MNOPO] array containing the current coordinates of the nodes.
	for (int i = 0; i != Geometry.getnumnodes(); i++) { this->file << std::scientific << Geometry.getintnode(i).getdx() << " " << std::scientific << Geometry.getintnode(i).getdy() << "\n"; }

	this->file << "/YD/YDN/D2NCI "		<< 21 << " " << Geometry.getdim() << " " << Geometry.getnumnodes()	<< "\n"; //[MNODIM][MNOPO] array containing the initial coordinates of the nodes(usually corresponding to undeformed system)
	for (int i = 0; i != Geometry.getnumnodes(); i++) { this->file << std::scientific << Geometry.getintnode(i).getdx() << " " << std::scientific << Geometry.getintnode(i).getdy() << "\n"; }
	
	this->file << "/YD/YDN/D2NFC "		<< 21 << " " << Geometry.getdim() << " " << Geometry.getnumnodes()	<< "\n"; //[MNODIM][MNOPO] array containing the current nodal forces.
	for (int i = 0; i != Geometry.getnumnodes(); i++) { this->file << std::scientific << Geometry.getintnode(i).getfx() << " " << std::scientific << Geometry.getintnode(i).getfy() << "\n"; }
	
	this->file << "/YD/YDN/D1NMCT 0\n"; //[MNOPO] array containing the current mass for each node.
	
	this->file << "/YD/YDN/D2NVC "		<< Geometry.getdim() << " " << Geometry.getnumnodes()	<< "\n";; //[MNODIM][MNOPO] array containing the current velocities for each node.
	for (int i = 0; i != Geometry.getnumnodes(); i++) { this->file << std::scientific << Geometry.getintnode(i).getvx() << " " << std::scientific << Geometry.getintnode(i).getvy() << "\n"; }
	
	this->file << "/YD/YDN/I1NOBF " << Geometry.getnumnodes() << "\n"; //[MNOPO] array containing a flag indicating that a node is on the boundary (usually set to 1 for all nodes indicating the boundary. Any node with the flag set to zero is considered not to be on the boundary)
	for (int i = 0; i != Geometry.getnumnodes(); i++) { this->file << Geometry.getintnode(i).getboundaryflag() << "\n"; }

	this->file << "/YD/YDN/I1NOPR " << Geometry.getnumnodes() << "\n"; //[MNOPO] array containing an ID of a property set associated with each node.
	for (int i = 0; i != Geometry.getnumnodes(); i++) { this->file << Geometry.getintnode(i).getmatid() << "\n"; }

	this->file.close();
	}
	else { std::cerr << "Error (in writenodes): Unable to open " << this->getyname() << std::endl; }
}

void Output::writeproperties(PropertyVariables &Property)
{
	this->file.open(this->getyname(), std::ios::out | std::ios::app | std::ios::binary); ;
	if (this->file.is_open())
	{
		this->file << "/YD/YDP/MPROP " << Property.getmaxnumprop() << "\n"; //Maximum number of properties sets.
		this->file << "/YD/YDP/NPROP " << Property.getnumprop() << "\n"; //Actual number of properties sets.

		this->file << "/YD/YDP/D1PCOH " << Property.getnumprop() << "\n"; //[MPROP] array containing internal cohesion
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getinternalcohesion() << " "; }

		this->file << "\n/YD/YDP/D1PICF " << Property.getnumprop() << "\n";//[MPROP] array containing interface friction.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getinterfacefriction() << " "; }
		
		this->file << "\n/YD/YDP/D1PEFR " << Property.getnumprop() << "\n";//[MPROP] array containing sliding friction (from 0 – glass to 2.5 – rubber).
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getslidingfriction() << " "; }
		
		this->file << "\n/YD/YDP/D1PEFT " << Property.getnumprop() << "\n";//[MPROP] array containing tensile strength of material for finite elements for each of the properties.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).gettensilestrength() << " "; }

		this->file << "\n /YD/YDP/D1PEGT " << Property.getnumprop() << "\n"; //[MPROP] array containing mode 1 energy rate
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getmodeienergyrate() << " "; }
		
		this->file << "\n/YD/YDP/D1PEGS " << Property.getnumprop() << "\n"; //[MPROP] array containing damping coefficient
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getmassdamping() << " "; }

		this->file << "\n/YD/YDP/D1PEKS " << Property.getnumprop() << "\n"; //[MPROP] array containing viscous damping of material for finite elements for each of the properties(for D1PEKS equal 2h E rho, finite element smaller than h is critically damped)
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getviscousdamping() << " "; }

		this->file << "\n/YD/YDP/D1PELA " << Property.getnumprop() << "\n";//[MPROP] array containing Lamé elastic constant. (lambda)
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getlambda() << " "; }

		this->file << "\n/YD/YDP/D1PEMU " << Property.getnumprop() << "\n";//[MPROP] array containing Lamé elastic constant. (mu)
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getmu() << " "; }

		this->file << "\n/YD/YDP/D1PEPE " << Property.getnumprop() << "\n";//[MPROP] array containing elastic penalty term for contact (usually 2 to 100 times greater than lambda)
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getelasticpenalty() << " "; }

		this->file << "\n/YD/YDP/D1PEPC " << Property.getnumprop() << "\n";//[MPROP] array containing contact penalty term for contact (usually 2 to 100 times greater than lambda)
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getcontactpenalty() << " "; }

		this->file << "\n/YD/YDP/D1PEPF " << Property.getnumprop() << "\n"; //[MPROP] array containing static friction coefficient 
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getinternalfriction() << " "; }

		this->file << "\n/YD/YDP/D1PBIF " << Property.getnumprop() << "\n"; //[MPROP] array containing interface friction coefficient
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getinterfacefriction() << " "; }

		this->file << "\n/YD/YDP/D1PERO " << Property.getnumprop() << "\n";//[MPROP] array containing density.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getdensity() << " "; }

		this->file << "\n/YD/YDP/D1PNAP " << Property.getnumprop() << "\n";//[MPROP] array containing amplitude of load applied as element surface pressure.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getloadamplitude() << " "; }

		this->file << "\n/YD/YDP/D1PNAF " << Property.getnumprop() << "\n";//[MPROP] array containing amplitude factor – all amplitudes are multiplied by this factor.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getamplitude() << " "; }

		this->file << "\n/YD/YDP/D1PNAI " << Property.getnumprop() << "\n";//[MPROP] array containing the increment of the amplitude factor each time step 
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getamplitudeinc() << " "; }

		this->file << "\n/YD/YDP/D1PNAT " << Property.getnumprop() << "\n";//[MPROP] array containing amplitude of load applied as element surface traction.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getloadamplitudeinc() << " "; }

		this->file << "\n/YD/YDP/D1PNAX " << Property.getnumprop() << "\n";//[MPROP] array containing amplitude of force or velocity in local x direction for each node.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getforcevelamplitudex() << " "; }

		this->file << "\n/YD/YDP/D1PNAY " << Property.getnumprop() << "\n";//[MPROP] array containing amplitude of force or velocity in local y direction for each node.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getforcevelamplitudey() << " "; }

		this->file << "\n/YD/YDP/D1PNXX " << Property.getnumprop() << "\n";//[MPROP] array containing the x component of x axis of local coordinate system for each node.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << 1 << " "; }

		this->file << "\n/YD/YDP/D1PNXY " << Property.getnumprop() << "\n";//[MPROP] array containing the y component of x axis of local coordinate system for each node.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << 0 << " "; }

		this->file << "\n/YD/YDP/D1PNYX " << Property.getnumprop() << "\n";//[MPROP] array containing the x component of y axis of local coordinate system for each node.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << 0 << " "; }

		this->file << "\n/YD/YDP/D1PNYY " << Property.getnumprop() << "\n";//[MPROP] array containing the y component of y axis of local coordinate system for each node.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << 1 << " "; }

		this->file << "\n/YD/YDP/D1PSEM " << Property.getnumprop() << "\n";//[MPROP] array containing the maximum tensile stretch.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << 0 << " "; }

		this->file << "\n/YD/YDP/D1PJRC " << Property.getnumprop() << "\n";//[MPROP] array containing joint JRC
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getjrc0() << " "; }

		this->file << "\n/YD/YDP/D1PJCS " << Property.getnumprop() << "\n";//[MPROP] array containing joint JCS
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getjcs0() << " "; }

		this->file << "\n/YD/YDP/D1PJSL " << Property.getnumprop() << "\n";//[MPROP] array containing joint length
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getjointsamplesize() << " "; }

		this->file << "\n/YD/YDP/I1PEFR " << Property.getnumprop() << "\n";//[MPROP] array containing fracture flag. If greater than zero all elements corresponding to this property will be processed as fracturing medium. If it is set to zero for a particular property set, no fracture processing for that property set will take place, i.e. finite elements that are assigned that property will remain non fractured.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getfractureflag() << " "; }

		this->file << "\n/YD/YDP/I1PEJP " << Property.getnumprop() << "\n";//[MPROP] array containing joint property(if > 0 joint)
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getjointproperty() << " "; }

		this->file << "\n/YD/YDP/I1PEMB " << Property.getnumprop() << "\n";//[MPROP] array. If greater than zero boundary nodes are marked, otherwise not.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << 0 << " "; }

		this->file << "\n/YD/YDP/I1PEMN " << Property.getnumprop() << "\n";//Numbers of successive mesh refinements.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << 0 << " "; }

		this->file << "\n/YD/YDP/I1PNFX " << Property.getnumprop() << "\n";//[MPROP] array. 1 indicates that a force is applied to a node in local x direction, 2 indicates that acceleration is applied to a node and 3 indicates that a velocity is applied to a node. For each node actual force, velocity or acceleration is obtained by multiplying amplitude (such as D1PNAX) with factor (such as D1PNAF).
		for (int i = 0; i != Property.getnummats(); i++) { this->file << 0 << " "; }

		this->file << "\n/YD/YDP/I1PNFY " << Property.getnumprop() << "\n";//[MPROP] array. 1 indicates that a force is applied to a node in local y direction, 2 indicates that acceleration is applied to a node and 3 indicates that a velocity is applied to a node. For each node actual force, velocity or acceleration is obtained by multiplying amplitude (such as D1PNAY) with factor (such as D1PNAF).
		for (int i = 0; i != Property.getnummats(); i++) { this->file << 0 << " "; }

		this->file << "\n/YD/YDP/I1PSDE " << Property.getnumprop() << "\n";//[MPROP] array containing ID of elastic damage state variable, say zero will indicate that first state variable is elastic damaged, one will indicate that elastic damaged is second state variable, i.e. stored in the second column of the array D2PELST.
		for (int i = 0; i != Property.getnummats(); i++) { this->file << 0 << " "; }

		this->file << "\n/YD/YDP/I1PTYP " << Property.getnumprop() << "\n";//[MPROP] array containing the type of each property. Indicates type of object to which this property is associated.For instance : (-1) 2D mechanical x, y d.o.f.node; (1) plane stress elastic triangle; (2) rigid triangle; (3) joint; (4) plain stress softening triangle
		for (int i = 0; i != Property.getnummats(); i++) { this->file << Property.getintmat(i).getpropertytype() << " "; }
		this->file << "\n";

		this->file.close();
	}
	else { std::cerr << "Error (in writeproperties): Unable to open " << this->getyname() << std::endl; }
}

void Output::writeoutput(PropertyVariables &Property, GeometryVariables &Geometry, Element &E, Node &N, ProblemVariables &Problem)
{
	this->file.open(this->getyname(), std::ios::out | std::ios::app | std::ios::binary);
	if (this->file.is_open())
	{
		this->file << "/YD/YDO/MOHYS     1\n"; //Maximum number of history variables.
		this->file << "/YD/YDO/NOHYS     1\n"; //Actual number of history variables.
		this->file << "/YD/YDO/DOHYP +5.000000000000e-004\n";//Output history accuracy.
		this->file << "/YD/YDO/D1OHYC     1\n"; //[MOHYS] output history factor to scale time.
		this->file << "+1.000000000000e+000\n";
		this->file << "/YD/YDO/D1OHYF     1\n"; //[MOHYS] output history factor to scale state.
		this->file << "+1.000000000000e+000\n";
		this->file << "/YD/YDO/D1OHYS     1\n"; //[MOHYS] output history state.
		this->file << "+0.000000000000e+000\n";
		this->file << "/YD/YDO/D1OHYT     1\n"; //[MOHYS] output history time.
		this->file << "+0.000000000000e+000\n";
		this->file << "/YD/YDO/D1OHYX     1\n"; //[MOHYS] output history x coordinate of the point.
		this->file << "+0.000000000000e+000\n";
		this->file << "/YD/YDO/D1OHYY     1\n"; //[MOHYS] output history y coordinate of the point.
		this->file << "+0.000000000000e+000\n";
		this->file << "/YD/YDO/I1OHYT     1\n"; //[MOHYS] output history type, i.e. which variable:0 – default (nothing)
			//1 – stress sigmaxx		11 – isotropic elastic damage
			//2 – stress sigmaxy		12 – borehole mass
			//3 – stress sigmayy		13 – borehole pressure
			//7 – velocity				14 – borehole spec.volume
			//8 – velocity x			15 – total kinetic energy
			//9 – velocity y			16 – G2 pressure
			//10 – velocity z			17 – isotropic elastic damage
		this->file << "15\n";

		////Gen_out_only
		//this->file << "Gen_out_only " << Problem.getgenoutonly() << "\n";

		////Material_List
		//for (int i = 0; i != Property.getnummats(); i++)
		//{
		//	this->file << "MATERIAL_LIST " << Property.getintmat(i).getnumber() << "\n";
		//	this->file << std::scientific 
		//		<< std::scientific << Property.getintmat(i).getdensity() << " "
		//		<< std::scientific << Property.getintmat(i).getyoungsmod() << " "
		//		<< std::scientific << Property.getintmat(i).getnu() << " "
		//		<< std::scientific << Property.getintmat(i).getmassdamping() << " "
		//		<< std::scientific << Property.getintmat(i).getelasticpenalty() << " "
		//		<< std::scientific << Property.getintmat(i).getcontactpenalty() << " "
		//		<< std::scientific << Property.getintmat(i).getmodeienergyrate() << " "
		//		<< std::scientific << Property.getintmat(i).getmodeiienergyrate() << " "
		//		<< std::scientific << Property.getintmat(i).gettensilestrength() << " "
		//		<< std::scientific << Property.getintmat(i).getinternalfriction() << " "
		//		<< std::scientific << Property.getintmat(i).getinternalcohesion() << " "
		//		<< std::scientific << Property.getintmat(i).getporefluidpressure() << " "
		//		<< std::scientific << Property.getintmat(i).getjointfriction() << " "
		//		<< std::scientific << Property.getintmat(i).getjrc0() << " "
		//		<< std::scientific << Property.getintmat(i).getjcs0() << " "
		//		<< std::scientific << Property.getintmat(i).getjointsamplesize() << " "
		//		<< std::scientific << Property.getintmat(i).getinterfacefriction() << " "
		//		<< std::scientific << Property.getintmat(i).getproblems2d() << "\n";
		//}

		//this->file << "Element_List " << Geometry.getnumele() << "\n";
		//for (int i = 0; i != Geometry.getnumele(); i++)
		//{
		//	this->file << i << " "
		//		<< Geometry.getintele(i).getfracturetype() << " 0\n";
		//}
		//
		//this->file << "\nInitial_Data " << Geometry.getnuminit() << "\n";
		//for (int i = 0; i != Geometry.getnumnodes(); i++)
		//{
		//	if (Geometry.getintnode(i).getinitflag() == 1)
		//	{
		//		this->file << Geometry.getintnode(i).getnumber() << " "
		//			<< 1 << " "
		//			<< std::scientific << Geometry.getintnode(i).getvx() << " "
		//			<< std::scientific << Geometry.getintnode(i).getvy() << " "
		//			<< std::scientific << Geometry.getintnode(i).getfx() << " "
		//			<< std::scientific << Geometry.getintnode(i).getfy() << " "
		//			<< std::scientific << Geometry.getintnode(i).getdx() << " "
		//			<< std::scientific << Geometry.getintnode(i).getdy() << "\n";
		//	}
		//}
		
		//this->file << "\nConstraints_List " << Geometry.getnumconstraints() << "\n";
		//for (int i = 0; i != Geometry.getnumnodes(); i++)
		//{
		//	if (Geometry.getintnode(i).getcflag() == 1)
		//	{
		//		this->file << i << " "
		//			<< Geometry.getintnode(i).getfcx() << " "
		//			<< Geometry.getintnode(i).getfcy() << " "
		//			<< Geometry.getintnode(i).getdcx() << " "
		//			<< Geometry.getintnode(i).getdcy() << " "
		//			<< Geometry.getintnode(i).getdcx() << " Velocity "
		//			<< Geometry.getintnode(i).getvcx() << " Velocity "
		//			<< Geometry.getintnode(i).getvcy() << "\n";
		//	}
		//}

		this->file << "$YDOIT\n"; //It is a command to execute the program.
		this->file << "$YSTOP\n"; //It is a command to stop the program.
		this->file.close();
	}
	else { std::cerr << "Error (in writeoutput): Unable to open " << this->getyname() << std::endl; }
}

int main(int argc, char **argv)
{
	//str gmshpath = argv[1];

	//arugments from the cmd line: argv[0]=./Ygen, argv[1]=...
	//ready all input files

	Input In;
	Output Out;
	UserVariables User;
	Surface S;
	Element E;
	Node N;
	GeometryVariables Geometry;
	Material Mat;
	PropertyVariables Property;
	ProblemVariables Problem;
	Customise Custom;

	Custom.setnodefile("../user/internal/nodes.txt");
	Custom.setnodevfile("../user/internal/nodalvelocities.txt");
	Custom.setnodeffile("../user/internal/nodalforces.txt");
	Custom.setndcfile("../user/internal/displ_constraints.txt");
	Custom.setnvcfile("../user/internal/vel_constraints.txt");
	Custom.setnfcfile("../user/internal/force_constraints.txt");
	Custom.setinitfile("../user/internal/initialdata.txt");
	Custom.setelementfile("../user/internal/elements.txt");

	Geometry.setdim(2);								//2D
	Geometry.setmaxdim(Geometry.getdim() + 1);
	Geometry.setmaxstatevar(Geometry.getdim());
	Geometry.setstatevar(Geometry.getdim());

	//Read in user data
	In.readcustom(Custom);
	In.readuser(User, Custom);
	In.readgeometry(Geometry, S, N, Custom);

	//write .geo file
	str basename = User.getoutfile().erase(User.getoutfile().size() - 3);
	Out.setgeoname(User.getlocaldir() + "/user/internal" + basename + ".geo");
	Out.setmeshname(User.getlocaldir() + "/user/internal" + basename + ".msh");
	Out.writegeo(argc, argv);

	//run gmsh
	system(("gmsh " + Out.getgeoname() + " -2").c_str());

	//read in geometry from .geo file;
	In.readnodes(Geometry, N, Custom, Property);
	In.readnodalvelocities(Geometry, N, Custom);
	In.readnodalforces(Geometry, N, Custom);

	In.readdisplconstraints(Geometry, N, Custom);
	In.readforceconstraints(Geometry, N, Custom);
	In.readvelconstraints(Geometry, N, Custom);

	In.readelements(Geometry, E, Custom);

	//read in material
	In.readmaterial(Mat, Property, Problem, Custom, Geometry);

	//read in simulation control parameters (problem data)
	In.readproblem(Problem, Property, Custom, Geometry);

	//write the y file

	Out.setyname(User.getlocaldir() + "/gid" + basename + "/" + basename + ".y");
	Out.writecontrol(Problem);
	Out.writeelements(Geometry);
	Out.writeinteractions(Problem);
	Out.writenodes(Geometry);
	Out.writeproperties(Property);
	Out.writeoutput(Property, Geometry, E, N, Problem);

	system("pause");

	return 0;
}