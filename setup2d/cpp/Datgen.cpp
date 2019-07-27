#include <algorithm>

#include <iostream>
#include <iomanip>
#include <iterator>
#include <map>
#include <fstream>
#include <functional>

#include <string>
#include <vector>

typedef std::string str;
typedef std::map<str, Material> strmatmap;			//Maps Strings to Materials
typedef std::map<int, Material> intmatmap;			//Maps Ints to Materials
typedef std::map<int, str> intstrmap;				//Maps Ints to Strings (Material Names)
typedef std::multimap<int, int> intintmultimap;		//Maps Node Numbers to Element Numbers 
typedef std::vector<str> strvec;
typedef std::vector<int> intvec;
typedef std::vector<double> doublevec;
typedef intintmultimap::iterator intintmultimapit;
typedef std::map<int, Element> intelemap;
typedef std::map<int, Node> intnodemap;
typedef std::pair<int, int> intintpair;
typedef std::pair<int, intvec> intintvecpair;
typedef std::pair<int, Node> intnodepair;
typedef std::pair<int, Element> intelepair;

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

public:
	PropertyVariables() {};   // This is the constructor declaration
	~PropertyVariables() {};  // This is the destructor: declaration

	int getmaxnumprop() { return this->maxnumprop; }
	int getnumprop() { return this->numprop; }
	Material getstrmat(str name);
	Material getintmat(int number);
	str getmatname(int number);

	int getnummats() { return this->nummats; }
	double getdtcrit();
	doublevec getdtcrits();

	void setmaxnumprop(int maxnumprop) { this->maxnumprop = maxnumprop; }
	void setnumprop(int numprop) { this->numprop = numprop; }

	void setstrmat(str name, Material Mat) { this->strmat[name] = Mat; }
	void setintmat(int number, Material Mat) { this->intmat[number] = Mat; }
	void setmatname(int number, str name) { this->matname[number] = name; }

	void setnummats(int nummats) { this->nummats = nummats; }
	void setdtcrit(double dtcrit) { this->dtcrit = dtcrit; }
	void setdtcrits(doublevec dtcrits) { this->dtcrits = dtcrits; }

};

str PropertyVariables::getmatname(int number)
{
	for (intstrmap::iterator it = this->matname.begin(); it != this->matname.end(); it++)
	{
		if (it->first == number)
		{
			return it->second;
		}
	}

	std::cerr << "Error: String index out of bound. Name not found." << std::endl;
	return " ";
}

Material PropertyVariables::getintmat(int number)
{
	intmatmap::iterator it = this->intmat.begin();

	for (intmatmap::iterator it = this->intmat.begin(); it != this->intmat.end(); it++)
	{
		if (it->first == number)
		{
			return it->second;
		}
	}

	std::cerr << "Error: Material index out of bound. Number: " << number << std::endl;
	Material Error;
	return Error;
}

Material PropertyVariables::getstrmat(str name)
{
	for (strmatmap::iterator it = this->strmat.begin(); it!= this->strmat.end(); it++)
	{
		if (it->first == name)
		{
			return it->second;
		}
	}

	std::cerr << "Error: Material index out of bound." << std::endl;
	Material Error;
	return Error;
}

doublevec PropertyVariables::getdtcrits()
{
	this->dtcrits.clear();
	for (strmatmap::iterator it = this->strmat.begin(); it != this->strmat.end(); ++it) 
	{ 
		this->dtcrits.push_back(it->second.getdtcrit()); 
	}
	return this->dtcrits;
}


double PropertyVariables::getdtcrit()
{
	this->dtcrits = this->getdtcrits();
	this->dtcrit = this->dtcrits[std::distance(this->dtcrits.begin(), std::min_element(this->dtcrits.begin(), this->dtcrits.end()))];
	return this->dtcrit;
}

class BonusVariables : public Variables
{

private:

	strvec data;

	double minelev;
	double minedge;
	double simtime;
	int numoutfiles;
	double projectilemeshsize;
	double plymeshsize;
	double interlayermeshsize;

public:
	BonusVariables() {};   // This is the constructor declaration
	~BonusVariables() {};  // This is the destructor: declaration

	strvec getdata() { return this->data; }

	double getminelev() { return this->minelev; }
	double getminedge() { return this->minedge; }
	double getsimtime() { return this->simtime; }
	int getnumoutfiles() { return this->numoutfiles; }
	double getprojectilemeshsize() { return this->projectilemeshsize; }
	double getplymeshsize() { return this->plymeshsize; };
	double getinterlayermeshsize() { return this->interlayermeshsize; }

	void setdata(strvec data) { this->data = data; }

	void setminelev(double minelev) { this->minelev = minelev; }
	void setminedge(double minedge) { this->minedge = minedge; }
	void setsimtime(double simtime) { this->simtime = simtime; }
	void setnumoutfiles(int numoutfiles) { this->numoutfiles = numoutfiles; }
	void setprojectilemeshsize(double projectilemeshsize) { this->projectilemeshsize = projectilemeshsize; }
	void setplymeshsize(double plymeshsize) { this->plymeshsize = plymeshsize; }
	void setinterlayermeshsize(double interlayermeshsize) { this->interlayermeshsize = interlayermeshsize; }
};

class ProblemVariables: public Variables
{

private:

	strvec problemdata;

	int maxndt;
	double gravity;
	double dtcrit;
	double dt;
	int outfreq;
	int restartfreq;
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

	int currentndt;
	double coordsize;
	double currentt;
	int contcouples;

public:
	ProblemVariables() {};   // This is the constructor declaration
	~ProblemVariables() {};  // This is the destructor: declaration

	strvec getproblemdata() { return this->problemdata; }

	str getaccuracy() { return this->accuracy; }
	str getfricmodel() { return this->fricmodel; }
	str getinitialaperturecorr() { return this->initialaperturecorr; }

	int getmaxndt() { return this->maxndt; }
	double getgravity() { return this->gravity; };
	double getdt() { return this->dt; };
	double getdtcrit() { return this->dtcrit; };
	int getoutfreq() { return this->outfreq; };
	int getrestartfreq() { return this->restartfreq; };
	int getgravitysettlingstage() { return this->gravitysettlingstage; };
	int getloadrampingstage() { return this->loadrampingstage; };
	int getmaxdim() { return this->maxdim; };
	double getmaxforce() { return this->maxforce; };
	int getmaxv() { return this->maxv; };
	double getmaxstress() { return this->maxstress; };
	int getmaxdispl() { return this->maxdispl; };
	double getmaxjointaperture() { return this->maxjointaperture; };
	int getmaxcontcouples() { return this->maxcontcouples; };
	double getbuffersize() { return this->buffersize; };

	double getcurrentdt() { return this->currentndt; };
	double getcoordsize() { return this->coordsize; };
	double getcurrentt() { return this->currentt; };
	int getcontcouples() { return this->contcouples; };

	void setproblemdata(strvec problemdata) { this->problemdata = problemdata; }

	void setmaxndt(int maxndt) { this->maxndt = maxndt; }
	void setgravity(double gravity) { this->gravity = gravity; }
	void setdt(double dt) { this->dt = dt; }
	void setdtcrit(double dtcrit) { this->dtcrit = dtcrit; }
	void setoutfreq(int outfreq) { this->outfreq = outfreq; }
	void setrestartfreq(int restartfreq) { this->restartfreq = restartfreq; }
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

class Node : public GeometryVariables
{
private:
	int number;
	int dx;
	int dy;
	int vx;
	int vy;
	int fx;
	int fy;
	int matid;
	double fcx; //force constraint in x
	double fcy;
	double vcx; //velocity constraint
	double vcy;
	double dcx; //displacement constraint
	double dcy;
	int cflag; //constraintflag
	int initflag; //initial data flag;

public:
	int getnumber() { return this->number; }

	int getdx() { return this->dx; }
	int getdy() { return this->dy; }

	int getvx() { return this->vx; }
	int getvy() { return this->vy; }

	int getfx() { return this->fx; }
	int getfy() { return this->fy; }

	int getmatid() { return this->matid; }
	int getcflag() { return this->cflag; }

	int getnumber() { return this->number; }

	int getdcx() { return this->dcx; }
	int getdcy() { return this->dcy; }

	int getfcx() { return this->fcx; }
	int getfcy() { return this->fcy; }

	int getvcx() { return this->vcx; }
	int getvcy() { return this->vcy; }

	int getinitflag() { return this->initflag; }

	///////////////////////////////////
	///////////////////////////////////

	void setdx(int dx) { this->dx = dx; }
	void setdy(int dy) { this->dy = dy; }

	void setvx(int vx) { this->vx = vx; }
	void setvy(int vy) { this->vy = vy; }

	void setfx(int fx) { this->fx = fx; }
	void setfy(int fy) { this->fy = fy; }

	void setmatid(int matid) { this->matid = matid; }
	void setcflag(int cflag) { this->cflag = cflag; }

	void setnumber(int number) { this->number = number; }

	void setdcx(int dcx) { this->dcx = dcx; }
	void setdcy(int dcy) { this->dcy = dcy; }

	void setfcx(int fcx) { this->fcx = fcx; }
	void setfcy(int fcy) { this->fcy = fcy; }

	void setvcx(int vcx) { this->vcx = vcx; }
	void setvcy(int vcy) { this->vcy = vcy; }

	void setinitflag(int flag) { this->initflag = flag; }

};

class Element : public GeometryVariables
{
private:
	int number;
	int fractureflag;
	int fracturetype; //Fracture or Visco-Elastic

public:
	int getnumber() { return this->number; }
	int getfractureflag() { return this->fractureflag; }
	int getfracturetype() { return this->fracturetype; }

	void setnumber() { this->number = number; }
	void setfractureflag() { this->fractureflag = fractureflag; }
	void setfracturetype() { this->fracturetype = fracturetype; }
};

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

	strvec elements;
	strvec nodes;
	strvec initnodes;
	strvec initdata;
	strvec nforce;
	strvec nvel;

	intintmultimap nodeelemap; //maps nodes to elements
	intelemap elemap;
	intnodemap nodemap;

public:
	GeometryVariables() 
	{
		this->numinit = 0;
		this->numconstraints = 0;
	};

	~GeometryVariables() {};

	int getdim() { return this->dim; }
	int getmaxdim() { return this->maxdim; }

	strvec getelements() { return this->elements; }
	strvec getnodes() { return this->nodes; }
	strvec getinitnodes() { return this->initnodes; }
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

	intvec getnodeele(int nodenumber);
	Element getintele(int elenumber) 
	{
		if (elenumber <= this->getnumele() & elenumber > 0)
		{
			return this->elemap[elenumber];
		}
		else
		{
			std::cerr << "Error (in GeometryVariables): Element number out of bound.";
		}
	}
	Node getintnode(int nodenumber) 
	{ 
		if (nodenumber <= this->getnumnodes() & nodenumber > 0)
		{
			return this->nodemap[nodenumber];
		}
		else
		{
			std::cerr << "Error (in GeometryVariables): Node number out of bound.";
		} 
	}

	void setdim(int dim) { this->dim = dim; }
	void setmaxdim(int maxdim) { this->maxdim = maxdim; }

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
	void setelements(strvec elements) { this->elements = elements; }
	void setnodes(strvec nodes) { this->nodes = nodes; }
	void setinitnodes(strvec initnodes) { this->initnodes = initnodes; }
	void setnforce(strvec nforce) { this->nforce = nforce; }
	void setnvel(strvec nvel) { this->nvel = nvel; }

	void setnodeele(int nodenum, int elenum);
	void setintnode(int nodenumber, Node N) { this->nodemap[nodenumber] = N; }
	void setintele(int elenumber, Element E) { this->elemap[elenumber] = E; }
};

void GeometryVariables::setnodeele(int nodenum, int elenum)
{
	std::pair<intintmultimapit, intintmultimapit> elements = nodeelemap.equal_range(nodenum);
	int count = std::distance(elements.first, elements.second);

	if (count > this->getmaxnodesperfe())
	{
		std::cerr << "Error: Maximum number of nodes per FE exceeded.\n";
	}
	else
	{
		this->nodeelemap.insert(intintpair (nodenum, elenum));
	}
}

intvec GeometryVariables::getnodeele(int nodenum)
{
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

class Input
{
public:
	void readelements(GeometryVariables &Geometry, Element &E);
	void readnodes(GeometryVariables &Geometry, Node &N);
	void readuser(UserVariables &User);
	void readproblem(ProblemVariables &Problem, BonusVariables &Bonus, PropertyVariables &Property);
	void readbonus(BonusVariables &Bonus);
	void readmaterial(Material Mat, PropertyVariables &Property, BonusVariables &Bonus);
	void readnodalvelocities(GeometryVariables &Geometry, Node &N);
	void readnodalforces(GeometryVariables &Geometry, Node &N);
	void readdisplconstraints(GeometryVariables &Geometry, Node &N);
	void readforceconstraints(GeometryVariables &Geometry, Node &N);
	void readvelconstraints(GeometryVariables &Geometry, Node &N);
};

void Input::readdisplconstraints(GeometryVariables &Geometry, Node &N)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int colno(1);
	int numdisplconstraints(0);

	infile.open("displ_constraints.txt");
	if (infile.is_open())
	{
		while (getline(infile, line)) 
		{
			if (line.at(0) != '#' && !line.empty())
			{
				lines.push_back(line);
				numdisplconstraints++;

				str delimiter = " ";
				size_t pos = 0;
				str token;
				while ((pos = line.find(delimiter)) != str::npos)
				{
					token = line.substr(0, pos);
					switch (colno)
					{
					case(0): {N.setnumber(std::stoi(token)); break; }
					case(1):
					{
						N.setdcx(std::stod(token));
						Geometry.getintnode(N.getnumber()).setdcx(std::stoi(token)); 
						break; 
					}
					case(2): 
					{
						N.setdcy(std::stod(token));
						Geometry.getintnode(N.getnumber()).setdcy(std::stoi(token)); 
						break; 
					}
					}
					line.erase(0, pos + delimiter.length());
					colno++;
				}
				if (N.getdcx() != 0 || N.getdcy() != 0)
				{
					N.setcflag(1);
					Geometry.getintnode(N.getnumber()).setcflag(1);
				}
			}
		}
		Geometry.setnumconstraints(Geometry.getnumconstraints() + numdisplconstraints);
		infile.close();
	}
	else { std::cerr << "Unable to read displacement constraints" << std::endl; }
}

void Input::readvelconstraints(GeometryVariables &Geometry, Node &N)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int colno(1);
	int numvelconstraints(0);

	infile.open("vel_constraints.txt");
	if (infile.is_open())
	{
		while (getline(infile, line)) 
		{
			if (line.at(0)!='#' && !line.empty())
			{
				lines.push_back(line);
				numvelconstraints++;

				str delimiter = " ";
				size_t pos = 0;
				str token;
				while ((pos = line.find(delimiter)) != str::npos) 
				{
					token = line.substr(0, pos);
					switch (colno)
					{
					case(0): {N.setnumber(std::stoi(token)); break; }
					case(1):
					{
						N.setvcx(std::stod(token));
						Geometry.getintnode(N.getnumber()).setvcx(std::stod(token)); 
						break; 
					}
					case(2): 
					{
						N.setvcy(std::stod(token));
						Geometry.getintnode(N.getnumber()).setvcy(std::stod(token)); 
						break; 
					}
					}
					line.erase(0, pos + delimiter.length());
					colno++;
				}
				if (N.getvcx() != 0 || N.getvcy() != 0)
				{
					N.setcflag(1);
					Geometry.getintnode(N.getnumber()).setcflag(1);
				}
			}
		}
		Geometry.setnumconstraints(Geometry.getnumconstraints() + numvelconstraints);
		infile.close();
	}
	else { std::cerr << "Unable to read velocity constraints" << std::endl; }
}

void Input::readforceconstraints(GeometryVariables &Geometry, Node &N)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int colno(1);
	int numforceconstraints(0);

	infile.open("force_constraints.txt");
	if (infile.is_open())
	{
		while (getline(infile, line)) 
		{
			if (line.at(0) != '#' && !line.empty())
			{
				lines.push_back(line);
				numforceconstraints++;

				str delimiter = " ";
				size_t pos = 0;
				str token;
				while ((pos = line.find(delimiter)) != str::npos)
				{
					token = line.substr(0, pos);
					switch (colno)
					{
					case(0): {N.setnumber(std::stoi(token)); break; }
					case(1): 
					{ 
						N.setfcx(std::stod(token));
						Geometry.getintnode(N.getnumber()).setfcx(std::stod(token)); 
						break; 
					}
					case(2): 
					{ 
						N.setfcy(std::stod(token));
						Geometry.getintnode(N.getnumber()).setfcy(std::stod(token)); 
						break; 
					}
					}
					line.erase(0, pos + delimiter.length());
					colno++;
				}
				if (N.getfcx() != 0 || N.getfcy() != 0)
				{
					N.setcflag(1);
					Geometry.getintnode(N.getnumber()).setcflag(1);
				}
			}
		}
		Geometry.setnumconstraints(Geometry.getnumconstraints() + numforceconstraints);
		infile.close();
	}
	else { std::cerr << "Unable to read force constraints" << std::endl; }
}

void Input::readelements(GeometryVariables &Geometry, Element &E)
{
	str line;
	strvec lines;
	std::ifstream infile;

	infile.open("elements.txt");
	if (infile.is_open())
	{
		while (getline(infile, line)) {
			lines.push_back(line);
		}
		infile.close();
	}
	else { std::cerr << "Unable to read elements" << std::endl; }

	Geometry.setelements(lines);
	Geometry.setnumele(lines.size());
	Geometry.setmaxnumele(10 * Geometry.getnumele());
}

void Input::readnodes(GeometryVariables &Geometry, Node &N)
{
	str line;
	strvec lines;
	std::ifstream infile;

	infile.open("nodes.txt", std::ifstream::in);
	if (infile.is_open())
	{
		while (getline(infile, line)) {
			lines.push_back(line);
		}
		infile.close();

		Geometry.setnodes(lines);
		Geometry.setnumnodes(lines.size());
		Geometry.setmaxnumnodes(50 * Geometry.getnumnodes());
	}
	else { std::cerr << "Unable to read nodes" << std::endl; }
}

void Input::readnodalforces(GeometryVariables &Geometry, Node &N)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int colno(1);
	int numinitforces(0);

	infile.open("nodalforces.txt");
	if (infile.is_open())
	{
		while (getline(infile, line))
		{
			if (line.at(0) != '#' && !line.empty())
			{
				lines.push_back(line);
				numinitforces++;

				str delimiter = " ";
				size_t pos = 0;
				str token;
				while ((pos = line.find(delimiter)) != str::npos)
				{
					token = line.substr(0, pos);
					switch (colno)
					{
					case(0): {N.setnumber(std::stoi(token)); break; }
					case(1): { Geometry.getintnode(N.getnumber()).setfx(std::stod(token)); break; }
					case(2): {Geometry.getintnode(N.getnumber()).setfy(std::stod(token)); break; }
					}
					line.erase(0, pos + delimiter.length());
					colno++;
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
	else { std::cerr << "Unable to read force constraints" << std::endl; }
}

void Input::readnodalvelocities(GeometryVariables &Geometry, Node &N)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int colno(1);
	int numinitvel(0);

	infile.open("nodalvelocities.txt");
	if (infile.is_open())
	{
		while (getline(infile, line))
		{
			if (line.at(0) != '#' && !line.empty())
			{
				lines.push_back(line);
				numinitvel++;

				str delimiter = " ";
				size_t pos = 0;
				str token;
				while ((pos = line.find(delimiter)) != str::npos)
				{
					token = line.substr(0, pos);
					switch (colno)
					{
					case(0): {N.setnumber(std::stoi(token)); break; }
					case(1): { Geometry.getintnode(N.getnumber()).setvx(std::stod(token)); break; }
					case(2): {Geometry.getintnode(N.getnumber()).setvy(std::stod(token)); break; }
					}
					line.erase(0, pos + delimiter.length());
					colno++;
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
	else { std::cerr << "Unable to read force constraints" << std::endl; }
}

void Input::readuser(UserVariables &User)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int linenumber(0);

	infile.open(User.getinfile(), std::ifstream::in);
	if (infile.is_open()) 
	{
		while (getline(infile, line)) 
		{
			if (line.at(0) != '#' && !line.empty())
			{

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
				lines.push_back(line);
				linenumber++;
			}
		}

		User.setuserdata(lines);
		infile.close();
	}
	else { std::cerr << "In readuser: Unable to open config.txt" << std::endl; }
}

void Input::readproblem(ProblemVariables &Problem, BonusVariables &Bonus, PropertyVariables &Property)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int linenumber(0);

	infile.open("problemdata.txt", std::ifstream::in);
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
				switch (linenumber)
				{
				case(0): {Problem.setmaxndt(std::stoi(line)); break; }
				case(1): {Problem.setgravity(std::stod(line)); break; }
				case(2): {Problem.setoutfreq(std::stoi(line)); break; }
				case(3): {Problem.setrestartfreq(std::stoi(line)); break; }
				case(4): {Problem.setgravitysettlingstage(std::stoi(line)); break; }
				case(5): {Problem.setloadrampingstage(std::stoi(line)); break; }
				case(6): {Problem.setmaxdim(std::stoi(line)); break; }
				case(7): {Problem.setmaxforce(std::stod(line)); break; }
				case(8): {Problem.setmaxv(std::stoi(line)); break; }
				case(9): {Problem.setmaxstress(std::stod(line)); break; }
				case(10): {Problem.setmaxdispl(std::stoi(line)); break; }
				case(11): {Problem.setmaxjointaperture(std::stod(line)); break; }
				case(12): {Problem.setmaxcontcouples(std::stoi(line)); break; }
				case(13): {Problem.setaccuracy(line); break; }
				case(14): {Problem.setfricmodel(line); break; }
				case(15): {Problem.setinitialaperturecorr(line); break; }
				}

				lines.push_back(line);
				linenumber++;
			}
		}
		infile.close();
	}
	else { std::cerr << "In readproblem: Unable to read problem parametersun" << std::endl; }

	Problem.setproblemdata(lines);

	Problem.setdt(std::floor(Property.getdtcrit()));

	Problem.setcurrentndt(0);
	Problem.setcoordsize(1);
	Problem.setcurrentt((double)0);
	Problem.setcontcouples(0);

	Problem.setbuffersize(0.1 * Bonus.getminedge());
}

void Input::readbonus(BonusVariables &Bonus)
{
	str line;
	strvec lines;
	std::ifstream infile;
	int linenumber(0);
	str delimiter = "//";

	infile.open("bonusvariables.txt", std::ifstream::in);
	if (infile.is_open())
	{
		while (getline(infile, line))
		{

			if (line.empty())
			{
				getline(infile, line);
			}

			if (line.at(0) != '#')
			{
				switch (linenumber)
				{
				case(0): {Bonus.setminelev(std::stod(line)); break; }
				case(1): {Bonus.setminedge(std::stod(line)); break; }
				case(2): {Bonus.setnumoutfiles(std::stoi(line)); break; }
				case(3): {Bonus.setprojectilemeshsize(std::stod(line)); break; }
				case(4): {Bonus.setplymeshsize(std::stod(line)); break; }
				case(5): {Bonus.setinterlayermeshsize(std::stod(line)); break; }
				}
				lines.push_back(line);
				linenumber++;
			}
		}
		Bonus.setdata(lines);
		infile.close();
	}
	else { std::cerr << "Unable to read bonus parameters" << std::endl; }
};

void Input::readmaterial(Material Mat, PropertyVariables &Property, BonusVariables &Bonus)
{
	str line;
	int linenumber(0);
	std::ifstream infile;
	Mat.setnumber(0);
	bool resetmat = false;

	infile.open("mat.txt", std::ifstream::in);

	if (infile.is_open())
	{
		while (getline(infile, line))
		{

			if (line.empty() || infile.bad())
			{
				linenumber = 0;
				getline(infile, line);
			}

			if (line.at(0) != '#')
			{
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

					Mat.setnumber(Mat.getnumber() + 1);
					break;
				}

				case(1): { Mat.setdensity(std::stod(line)); break; }
				case(2): { Mat.setyoungsmod(std::stod(line)); break; }
				case(3): 
				{ 
					Mat.setnu(std::stod(line)); 
					Mat.setlambda(Mat.getnu() * Mat.getyoungsmod() / (1 + Mat.getnu()) / (1 - 2 * Mat.getnu()));
					Mat.setmu(Mat.getyoungsmod() / (2 * (1 + Mat.getnu()))); break; 
				}
				case(4): { Mat.setmassdamping(std::stod(line)); break; }
				case(5): 
				{ 
					Mat.setelasticpenalty(std::stod(line)); 
					Mat.setdtcrit(Bonus.getminelev() * Mat.getdensity() / Mat.getelasticpenalty()); 
					break;
				}
				case(6): { Mat.setcontactpenalty(std::stod(line)); break; }
				case(7): { Mat.setmodeienergyrate(std::stod(line)); break; }
				case(8): { Mat.setmodeiienergyrate(std::stod(line)); break; }
				case(9): { Mat.settensilestrength(std::stod(line)); break; }
				case(10): { Mat.setinternalfriction(std::stod(line)); break; }
				case(11): { Mat.setinternalcohesion(std::stod(line)); break; }
				case(12): { Mat.setporefluidpressure(std::stod(line)); break; }
				case(13): { Mat.setjointfriction(std::stod(line)); break; }
				case(14): { Mat.setjrc0(std::stod(line)); break; }
				case(15): { Mat.setjcs0(std::stod(line)); break; }
				case(16): { Mat.setjointsamplesize(std::stod(line)); break; }
				case(17): { Mat.setinterfacefriction(std::stod(line)); break; }
				case(18): { Mat.setproblems2d(line); break; }
				case(19): { Mat.setfractureflag(stoi(line)); break; }
				case(20): { Mat.setnumberofmeshrefinements(stoi(line)); break; }
				case(21): { Mat.setpropertytype(stoi(line)); break; }
				case(22): { Mat.setslidingfriction(stod(line)); break; }
				case(23): { Mat.setjointproperty(stoi(line)); break; }
				}
				
				linenumber++;
			}
		}

		infile.close();
	}
	else { std::cerr << "Unable to open mat.txt" << std::endl; }

	Property.setnummats(Mat.getnumber());
	Property.setnumprop(Property.getnummats());
	Property.setmaxnumprop(Property.getnummats()+1);
};

class Output
{
	friend class Material;
	friend class PropertyVariables;
public:

	Output() {};
	~Output() {};

	std::ofstream file;
	str name;

	str getname() { return this->name; }
	void setname(str name) { this->name = name; }

	void writevector(const strvec &lines);
	void writecontrol(ProblemVariables &Problem, BonusVariables &Bonus);
	void writeelements(GeometryVariables &Geometry);
	void writeinteractions(ProblemVariables &Problem);
	void writenodes(GeometryVariables &Geometry);
	void writeproperties(PropertyVariables &Property);
	void writeoutput(PropertyVariables &Property, GeometryVariables &Geometry, Element &Ele, Node &N);
};		

void Output::writevector(const strvec &lines)
{
	this->file.open(this->getname(), std::ios::out | std::ios::app | std::ios::binary); ;
	if (this->file.is_open())
	{
		if (!lines.empty())
		{
			//std::copy(lines.begin(), lines.end(), std::ostream_iterator<str>(this->file, " "));
			//this->file << lines.back() << "\n";
			for (strvec::const_iterator it = lines.begin(); it != lines.end(); it++)
			{
				this->file << *it << " ";
			}

		}
		this->file.close();
	}
	else { std::cerr << "In writevector: Unable to open output file." << std::endl; }
}

void Output::writecontrol(ProblemVariables &Problem, BonusVariables &Bonus)
{
	this->file.open(this->getname(), std::ios::out | std::ios::trunc | std::ios::binary); ;
	if (this->file.is_open())
	{
		this->file << "\t/YD/YDC/MCSTEP "		 				<< Problem.getmaxndt()					<< "\n"; //Maximum number of steps
		this->file << "\t\t/YD/YDC/NCSTEP "						<< Problem.getcurrentdt()				<< "\n"; //Current / Actual number of time steps (cannot be greater than /YD/YDC/MCSTEP)
		this->file << "\t\t/YD/YDC/DCGRAY "						<< Problem.getgravity()					<< "\n"; //Acceleration of gravity (in y direction)
		this->file << "\t\t/YD/YDC/DCSTEC " << std::scientific	<< Problem.getdt()						<< "\n"; //Size of the time step (see in the book how to calculate critical-maximum time step)
		this->file << "\t\t/YD/YDC/ICOUTF "						<< Problem.getoutfreq()					<< "\n"; //Output frequency – every so many time steps complete state of the system is recorded in a file with extension .ym which can be visualized using M program, which is FEM/DEM Visualizer accompanying Y program.
		this->file << "\t\t/YD/YDC/ICSAVF "						<< Problem.getrestartfreq()				<< "\n"; //Restart save frequency
		this->file << "\t\t/YD/YDC/DCGRST "						<< Problem.getgravitysettlingstage()	<< "\n"; //Gravity settling stage
		this->file << "\t\t/YD/YDC/DCRMPT "						<< Problem.getloadrampingstage()		<< "\n"; //Load ramping stage
		this->file << "\t\t/YD/YDC/DCSIZC "						<< Problem.getcoordsize()				<< "\n"; //Maximum size of coordinate in any direction (size of physical space – corresponding outputs are normalized using this value)
		this->file << "\t\t/YD/YDC/DCSIZF " << std::scientific	<< Problem.getmaxforce()				<< "\n"; //Maximum size of force in any direction(corresponding outputs are normalized using this value)
		this->file << "\t\t/YD/YDC/DCSIZV "						<< Problem.getmaxv()					<< "\n"; //Maximum size of velocity in any direction (corresponding outputs are normalized using this value)
		this->file << "\t\t/YD/YDC/DCSIZS " << std::scientific	<< Problem.getmaxstress()				<< "\n"; //Maximum stress (buffer size)
		this->file << "\t\t/YD/YDC/DCSIZD "						<< Problem.getmaxdispl()				<< "\n"; //Maximum size of displacement in any direction
		this->file << "\t\t/YD/YDC/DCSIZA " << std::scientific	<< Problem.getmaxjointaperture()		<< "\n"; //Maximum joint aperture
		this->file << "\t\t/YD/YDC/DCTIME "						<< Problem.getcurrentt()				<< "\n"; //Current time. i.e. time at start of this run.
		this->file << "\t\t/YD/YDC/ICOUTI "						<< 0									<< "\n"; //Current number of iterations; Serial number of first output associated with this run.
		this->file << "\t\t/YD/YDC/ICOUTP "						<< 4									<< "\n"; //Number of characters for each number in output file (for example: three characters is equivalent to six significant digits)
		this->file << "\t\t/YD/YDC/ICFMTY "						<< 0									<< "\n"; //
		this->file << "\t\t/YD/YDC/ICIATY "						<< 0									<< "\n\n"; //
		this->file.close();
	}
	else { std::cerr << "In writecontrol: Unable to open " << this->getname() << std::endl; }
}

void Output::writeelements(GeometryVariables &Geometry)
{
	this->file.open(this->getname(), std::ios::out | std::ios::app | std::ios::binary); ;
	if (this->file.is_open())
	{
		this->file << "/YD/YDE/MELEM "  << Geometry.getmaxnumele()									<< "\n"; //Maximum number of finite elements(with fracture the actual number of finite elements increases during the run, however, it should not exceed this number)
		this->file << "/YD/YDE/NELEM "  << Geometry.getnumele()										<< "\n"; //Actual number of finite elements at the beginning of this run.
		this->file << "/YD/YDE/MELST "  << Geometry.getmaxstatevar()								<< "\n"; //Maximum number of state variables per finite element.
		this->file << "/YD/YDE/NELST "  << Geometry.getstatevar()									<< "\n"; //Actual number of state variables per finite element.
		this->file << "/YD/YDE/MELNO "  << Geometry.getmaxnodesperfe()								<< "\n"; //Maximum number of nodes per finite element.
		this->file << "/YD/YDE/NELNO "  << Geometry.getnodesperfe()									<< "\n"; //Actual number of nodes per finite element.
		this->file << "/YD/YDE/D2ELST " << Geometry.getstatevar() << " " << Geometry.getnumele()	<< "\n"; //[MELST][MELEM] array containing state variables for all finite elements.
		this->file << "/YD/YDE/I1ELCF " << Geometry.getnumele()										<< "\n"; //Head of a list of contacting couples for each finite element.
		this->file << "/YD/YDE/I1ELPR " << Geometry.getnumele()										<< "\n"; //Set of properties associated with each element.
		for (int i = 0; i < Geometry.getnumele(); i++) 
		{ 
			this->file << "1\n";
		};
		this->file << "\n/YD/YDE/I2ELTO " << 21 << " " << 3 << " " << Geometry.getnumele()	<< "\n"; //[MELNO][MELEM] topology array containing nodes for each finite element.
		this->writevector(Geometry.getelements());
		//this->file << "/YD/YDE/I1ELTY " << Geometry.getnumele() << "\n";
		//for (int i = 0; i < Geometry.getnumele(); i++) { this->file << "-1 "; };
		//this->file << "\n";
		//this->file << "\t/YD/YDE/I2ELJP   21  3  0\n"; 
		this->file << "\n";
		this->file.close();
	}
	else { std::cerr << "In writeelements: Unable to open " << this->getname() << std::endl; }
}

void Output::writeinteractions(ProblemVariables &Problem)
{
	this->file.open(this->getname(), std::ios::out | std::ios::app | std::ios::binary); ;
	if (this->file.is_open())
	{
		this->file << "/YD/YDI/MICOUP "	<< Problem.getmaxcontcouples()	<< "\n"; //Maximum number of contacting couples of finite elements.
		this->file << "/YD/YDI/NICOUP "	<< Problem.getcontcouples()		<< "\n"; //Actual number of contacting couples of finite elements (always set to zero)
		this->file << "/YD/YDI/IIECFF "	<< -2							<< "\n"; //Internal variable used for contact (always set to -2)
		this->file << "/YD/YDI/DIEDI "	<< 200							<< "\n"; //Internal variable (travel since last detection) used to trigger contact detection (always set to a large number in order to trigger contact detection immediately at the start of this run)
		this->file << "/YD/YDI/DIEZON "	<< Problem.getbuffersize()		<< "\n"; //Contact detection buffer size, Size of the buffer around each finite element for contact detection purposes (usually 1/5 of the size of the smallest finite element)
		this->file << "/YD/YDI/D1IESL "	<< 0							<< "\n"; //[MELEM] array used to store sliding distance between couples in contact
		this->file << "/YD/YDI/I1IECN "	<< 0							<< "\n"; //[MICOUP] array used to store next contacting couple in the list.
		this->file << "/YD/YDI/I1IECT "	<< 0							<< "\n"; //[MICOUP] array used to store target finite element for each contacting couple.
		this->file << "\n";
		this->file.close();
	}
	else
	{ std::cerr << "In writeinteractions: Unable to open " << this->getname() << std::endl; }
}

void Output::writenodes(GeometryVariables &Geometry)
{
	this->file.open(this->getname(), std::ios::out | std::ios::app | std::ios::binary); ;
	if (this->file.is_open())
	{
	this->file << "/YD/YDN/MNODIM "		<< Geometry.getmaxdim()												<< "\n";	//Maximum number of degrees of freedom per node (usually 2 in 2D)
	this->file << "/YD/YDN/NNODIM "		<< Geometry.getdim()												<< "\n";	//Actual number of degrees of freedom per node (usually 2 in 2D)
	this->file << "/YD/YDN/MNOPO "		<< Geometry.getmaxnumnodes()										<< "\n";	//Maximum number of nodes.
	this->file << "/YD/YDN/NNOPO "		<< Geometry.getnumnodes()											<< "\n\n";	//Actual number of nodes.
	
	this->file << "/YD/YDN/D2NCC "		<< Geometry.getdim() << " " << Geometry.getnumnodes()	<< "\n"; //[MNODIM][MNOPO] array containing the current coordinates of the nodes.
	this->writevector(Geometry.getnodes());

	this->file << "/YD/YDN/D2NCI "		<< Geometry.getdim() << " " << Geometry.getnumnodes()	<< "\n"; //[MNODIM][MNOPO] array containing the initial coordinates of the nodes(usually corresponding to undeformed system)
	this->writevector(Geometry.getnodes());
	
	this->file << "/YD/YDN/D2NFC "		<< Geometry.getdim() << " " << Geometry.getnumnodes()	<< "\n"; //[MNODIM][MNOPO] array containing the current nodal forces.
	this->writevector(Geometry.getnforce());
	
	this->file << "/YD/YDN/D1NMCT "		<< 0 << "\n\n"; //[MNOPO] array containing the current mass for each node.
	
	this->file << "/YD/YDN/D2NVC "		<< 21 << " " << Geometry.getdim() << " " << Geometry.getnumnodes()	<< "\n";; //[MNODIM][MNOPO] array containing the current velocities for each node.
	this->writevector(Geometry.getnvel());
	
	this->file << "/YD/YDN/I1NOBF " << Geometry.getnumnodes() << " " << 0 << "\n"; //[MNOPO] array containing a flag indicating that a node is on the boundary (usually set to 1 for all nodes indicating the boundary. Any node with the flag set to zero is considered not to be on the boundary)
	
	this->file << "/YD/YDN/I1NOPR " << Geometry.getnumnodes() << " " << 0 << "\n"; //[MNOPO] array containing an ID of a property set associated with each node.
	
	this->file.close();
	}
	else { std::cerr << "In writenodes: Unable to open " << this->getname() << std::endl; }
}

void Output::writeproperties(PropertyVariables &Property)
{
	this->file.open(this->getname(), std::ios::out | std::ios::app | std::ios::binary); ;
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
	else { std::cerr << "In writeproperties: Unable to open " << this->getname() << std::endl; }
}

void Output::writeoutput(PropertyVariables &Property, GeometryVariables &Geometry, Element &E, Node &N)
{
	this->file.open(this->getname(), std::ios::out | std::ios::app | std::ios::binary);
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
		this->file << "15\n\n";

		//Gen_out_only
		this->file << "Gen_out_only Yes\n";

		//Material_List
		for (int i = 0; i != Property.getnummats(); i++)
		{
			this->file << "MATERIAL_LIST" << Property.getintmat(i).getnumber() << "\n";
			this->file << std::scientific 
				<< std::scientific << Property.getintmat(i).getdensity() << " "
				<< std::scientific << Property.getintmat(i).getyoungsmod() << " "
				<< std::scientific << Property.getintmat(i).getnu() << " "
				<< std::scientific << Property.getintmat(i).getmassdamping() << " "
				<< std::scientific << Property.getintmat(i).getelasticpenalty() << " "
				<< std::scientific << Property.getintmat(i).getcontactpenalty() << " "
				<< std::scientific << Property.getintmat(i).getmodeienergyrate() << " "
				<< std::scientific << Property.getintmat(i).getmodeiienergyrate() << " "
				<< std::scientific << Property.getintmat(i).gettensilestrength() << " "
				<< std::scientific << Property.getintmat(i).getinternalfriction() << " "
				<< std::scientific << Property.getintmat(i).getinternalcohesion() << " "
				<< std::scientific << Property.getintmat(i).getporefluidpressure() << " "
				<< std::scientific << Property.getintmat(i).getjointfriction() << " "
				<< std::scientific << Property.getintmat(i).getjrc0() << " "
				<< std::scientific << Property.getintmat(i).getjcs0() << " "
				<< std::scientific << Property.getintmat(i).getjointsamplesize() << " "
				<< std::scientific << Property.getintmat(i).getinterfacefriction() << " "
				<< std::scientific << Property.getintmat(i).getproblems2d() << " ";
		}

		this->file << "Element_List " << Geometry.getnumele() << "\n";
		for (int i = 0; i != Geometry.getnumele(); i++)
		{
			this->file << i << " "
				<< Geometry.getintele(i).getfracturetype() << " 0\n";
		}
		
		this->file << "\nInitial_Data " << Geometry.getnuminit() << "\n";
		for (int i = 0; i != Geometry.getnumnodes(); i++)
		{
			if (Geometry.getintnode(i).getinitflag() == 1)
			{
				this->file << Geometry.getintnode(i).getnumber() << " "
					<< 1 << " "
					<< std::scientific << Geometry.getintnode(i).getvx() << " "
					<< std::scientific << Geometry.getintnode(i).getvy() << " "
					<< std::scientific << Geometry.getintnode(i).getfx() << " "
					<< std::scientific << Geometry.getintnode(i).getfy() << " "
					<< std::scientific << Geometry.getintnode(i).getdx() << " "
					<< std::scientific << Geometry.getintnode(i).getdy() << "\n"
			}
		}
		
		this->file << "\nConstraints_List " << Geometry.getnumconstraints() << "\n";
		for (int i = 0; i != Geometry.getnumnodes(); i++)
		{
			if (Geometry.getintnode(i).getcflag() == 1)
			{
				this->file << i << " "
					<< Geometry.getintnode(i).getfcx() << " "
					<< Geometry.getintnode(i).getfcy << " "
					<< Geometry.getintnode(i).getdcx << " "
					<< Geometry.getintnode(i).getdcy << " "
					<< Geometry.getintnode(i).getdcx << " Velocity "
					<< Geometry.getintnode(i).getvcx << " Velocity "
					<< Geometry.getintnode(i).getvcy << "\n";
			}
		}

		this->file << "$YDOIT\n"; //It is a command to execute the program.
		this->file << "$YSTOP\n"; //It is a command to stop the program.
		this->file.close();
	}
	else { std::cerr << "In writeoutput: Unable to open " << this->getname() << std::endl; }
}

int main(int argc, char **argv)
{
	//arugments from the cmd line: argv[0]=./Ygen, argv[1]=...
	//ready all input files

	Input In;
	Output Out;
	UserVariables User;
	Element Ele;
	Node N;
	GeometryVariables Geometry;
	BonusVariables Bonus;
	Material Mat;
	PropertyVariables Property;
	ProblemVariables Problem;

	User.setinfile("config.txt");

	Geometry.setdim(2);								//2D
	Geometry.setmaxdim(Geometry.getdim() + 1);
	Geometry.setmaxstatevar(Geometry.getdim());
	Geometry.setstatevar(Geometry.getdim());

	In.readuser(User);
	In.readelements(Geometry, Ele);
	In.readnodes(Geometry, N);
	In.readdisplconstraints(Geometry, N);
	In.readforceconstraints(Geometry, N);
	In.readvelconstraints(Geometry, N);
	In.readbonus(Bonus);
	In.readmaterial(Mat, Property, Bonus);
	In.readproblem(Problem, Bonus, Property);

	Out.setname(User.getoutfile().erase(User.getoutfile().size() - 2) + ".y");

	Out.writecontrol(Problem, Bonus);
	Out.writeelements(Geometry);
	Out.writeinteractions(Problem);
	Out.writenodes(Geometry);
	Out.writeproperties(Property);
	Out.writeoutput(Property, Geometry, Ele, N);

	system("pause");

	return 0;
}