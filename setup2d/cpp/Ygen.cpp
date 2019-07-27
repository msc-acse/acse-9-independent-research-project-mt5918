#include <algorithm>

#include <iostream>
#include <iomanip>
#include <fstream>

#include <string>
#include <vector>

class Variables
{

};

class Material : public Variables
{

private:
	double density;
	double youngsmod;
	double poisson;
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
	std::string problems2d;

	std::vector<std::string> data;

public:
	Material(std::vector <std::string> lines);   // This is the constructor declaration
	~Material() {};  // This is the destructor: declaration

	double getdensity() { return this->density; }
	double getyoungsmod() { return this->youngsmod; };
	double getpoisson() { return this->poisson; };
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
	std::string getproblems2d() { return this->problems2d; };

	std::vector<std::string> getdata() { return this->data; };

	void setdensity(double density) { this->density = density; }
	void setyoungsmod(double youngsmod) { this->youngsmod = youngsmod; }
	void setpoisson(double poisson) { this->poisson = poisson; }
	void setmassdamping(double massdamping) { this->massdamping = massdamping; }
	void setelasticpenalty(double dtcritsteel) { this->elasticpenalty = elasticpenalty; }
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
	void setproblems2d(std::string problems2d) { this->problems2d = problems2d; }

	void setdata(std::vector<std::string> data) { this->data = data; }
};

Material::Material(std::vector <std::string> lines)
{
	Material::setdensity(std::stod(lines[19]));
	Material::setyoungsmod(std::stod(lines[20]));
	Material::setpoisson(std::stod(lines[21]));
	Material::setmassdamping(std::stod(lines[22]));
	Material::setelasticpenalty(std::stod(lines[23]));
	Material::setcontactpenalty(std::stod(lines[24]));
	Material::setmodeienergyrate(std::stod(lines[25]));
	Material::setmodeiienergyrate(std::stod(lines[26]));
	Material::settensilestrength(std::stod(lines[27]));
	Material::setinternalfriction(std::stod(lines[28]));
	Material::setinternalcohesion(std::stod(lines[29]));
	Material::setporefluidpressure(std::stod(lines[30]));
	Material::setjointfriction(std::stod(lines[31]));
	Material::setjrc0(std::stod(lines[32]));
	Material::setjcs0(std::stod(lines[33]));
	Material::setjointsamplesize(std::stod(lines[34]));
	Material::setinterfacefriction(std::stod(lines[35]));
	Material::setproblems2d(lines[36]);
}

class BonusVariables : public Variables
{

private:
	double minelev;
	double minedge;
	double simtime;
	double dtcritglas;
	double dtcritsteel;
	double dtcritpvb;
	double dtcritrock;
	int numoutfiles;
	double projectilemeshsize;
	double plymeshsize;
	double interlayermeshsize;

public:
	BonusVariables(std::vector <std::string> lines, Material Rock, Material Steel, Material Pvb, Material Glas);   // This is the constructor declaration
	~BonusVariables() {};  // This is the destructor: declaration

	double getminelev() { return this->minelev; }
	double getminedge() { return this->minedge; };
	double getsimtime() { return this->simtime; };
	double getdtcritglas() { return this->dtcritglas; };
	double getdtcritsteel() { return this->dtcritsteel; };
	double getdtcritpvb() { return this->dtcritpvb; };
	double getdtcritrock() { return this->dtcritrock; };
	int getnumoutfiles() { return this->numoutfiles; };
	double getprojectilemeshsize() { return this->projectilemeshsize; };
	double getplymeshsize() { return this->plymeshsize; };
	double getinterlayermeshsize() { return this->interlayermeshsize; };

	void setminelev(double maxndt) { this->minelev = minelev; }
	void setminedge(double minedge) { this->minedge = minedge; }
	void setsimtime(double simtime) { this->simtime = simtime; }
	void setdtcritglas(double dtcritglas) { this->dtcritglas = dtcritglas; }
	void setdtcritsteel(double dtcritsteel) { this->dtcritsteel = dtcritsteel; }
	void setdtcritpvb(double dtcritpvb) { this->dtcritpvb = dtcritpvb; }
	void setdtcritrock(double dtcritrock) { this->dtcritrock = dtcritrock; }
	void setnumoutfiles(int numoutfiles) { this->numoutfiles = numoutfiles; }
	void setprojectilemeshsize(double projectilemeshsize) { this->projectilemeshsize = projectilemeshsize; }
	void setplymeshsize(double plymeshsize) { this->plymeshsize = plymeshsize; }
	void setinterlayermeshsize(double interlayermeshsize) { this->interlayermeshsize = interlayermeshsize; }
};

BonusVariables::BonusVariables(std::vector <std::string> lines, Material Rock, Material Steel, Material Pvb, Material Glas)
{
	BonusVariables::setminelev(std::stod(lines[12]));
	BonusVariables::setminedge(std::stod(lines[13]));
	BonusVariables::setsimtime(std::stod(lines[14]));
	BonusVariables::setdtcritglas(std::stod(lines[15]));
	BonusVariables::setdtcritsteel(std::stod(lines[16]));
	BonusVariables::setdtcritpvb(std::stod(lines[17]));
	BonusVariables::setdtcritrock(std::stod(lines[18]));
	BonusVariables::setnumoutfiles(std::stoi(lines[19]));
	BonusVariables::setprojectilemeshsize(std::stod(lines[20]));
	BonusVariables::setplymeshsize(std::stod(lines[21]));
	BonusVariables::setinterlayermeshsize(std::stod(lines[22]));
}

class ProblemVariables: public Variables
{

private:
	int maxndt;
	double gravity;
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
	std::string accuracy;
	std::string fricmodel;
	std::string initialaperturecorr;

	int currentndt;
	double coordsize;
	double currentt;
	int contcouples;

public:
	ProblemVariables(std::vector <std::string> lines, BonusVariables Bonus);   // This is the constructor declaration
	~ProblemVariables() {};  // This is the destructor: declaration

	std::string getaccuracy() { return this->accuracy; }
	std::string getfricmodel() { return this->fricmodel; }
	std::string getinitialaperturecorr() { return this->initialaperturecorr; }

	int getmaxndt() { return this->maxndt; }
	double getgravity() { return this->gravity; };
	double getdt() { return this->dt; };
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

	void setmaxndt(int maxndt) { this->maxndt = maxndt; }
	void setgravity(double gravity) { this->gravity = gravity; }
	void setdt(double dt) { this->dt = dt; }
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
	void setaccuracy(std::string accuracy) { this->accuracy = accuracy; }
	void setfricmodel(std::string fricmodel) { this->fricmodel = fricmodel; }
	void setinitialaperturecorr(std::string initialaperturecorr) { this->initialaperturecorr = initialaperturecorr; }

	void setcurrentndt(int currentndt) { this->currentndt = currentndt; }
	void setcoordsize(int coordsize) { this->coordsize = coordsize; }
	void setcurrentt(double currentt) { this->currentt = currentt; }
	void setcontcouples(int contcouples) { this->contcouples = contcouples; }
};

ProblemVariables::ProblemVariables(std::vector <std::string> lines, BonusVariables Bonus)
{
	//Check and calculate problem data -this was already provisionally done in an excel sheet, now with code
	std::vector <double> dtcrits{ Bonus.getdtcritglas(), Bonus.getdtcritrock(), Bonus.getdtcritpvb(), Bonus.getdtcritsteel() };
	double dtcrit, dt, buffersize;

	std::vector<double>::iterator dtcrit = std::min_element(std::begin(dtcrits), std::end(dtcrits));
	dt = floor(dtcrit);
	buffersize = 0.1 * Bonus.getminedge();

	ProblemVariables::setmaxndt(std::stoi(lines[17]));
	ProblemVariables::setgravity(std::stod(lines[18]));
	ProblemVariables::setdt(dt);
	ProblemVariables::setoutfreq(std::stoi(lines[19]));
	ProblemVariables::setrestartfreq(std::stoi(lines[20]));
	ProblemVariables::setgravitysettlingstage(std::stoi(lines[21]));
	ProblemVariables::setloadrampingstage(std::stoi(lines[22]));
	ProblemVariables::setmaxdim(std::stoi(lines[23]));
	ProblemVariables::setmaxforce(std::stod(lines[24]));
	ProblemVariables::setmaxv(std::stoi(lines[25]));
	ProblemVariables::setmaxstress(std::stod(lines[26]));
	ProblemVariables::setmaxdispl(std::stoi(lines[27]));
	ProblemVariables::setmaxjointaperture(std::stod(lines[28]));
	ProblemVariables::setmaxcontcouples(std::stoi(lines[29]));
	ProblemVariables::setbuffersize(buffersize);
	ProblemVariables::setaccuracy(lines[30]);
	ProblemVariables::setfricmodel(lines[31]);
	ProblemVariables::setinitialaperturecorr(lines[33]);

	ProblemVariables::setcurrentndt(0);
	ProblemVariables::setcoordsize(1);
	ProblemVariables::setcurrentt((double)0);
	ProblemVariables::setcontcouples(0);
};

class UserVariables : public Variables
{
private:
	std::string makenewfile;
	std::string file;
	std::string login;
	std::string hpcdir;
	std::string localdir;
	std::string qsubdir;

public:
	UserVariables(std::vector <std::string> variables);   // This is the constructor declaration
	~UserVariables() {};  // This is the destructor: declaration

	std::string getmakenewfile() { return this->makenewfile; }
	std::string getfile() { return this->file; }
	std::string getlogin() { return this->login; }
	std::string gethpcdir() { return this->hpcdir; }
	std::string getlocaldir() { return this->localdir; }
	std::string getqsubdir() { return this->qsubdir; }

	void setmakenewfile(std::string makenewfile) { this->makenewfile = makenewfile; }
	void setfile(std::string file) { this->file = file; }
	void setlogin(std::string login) { this->login = login; }
	void sethpcdir(std::string hpcdir) { this->hpcdir = hpcdir; }
	void setlocaldir(std::string localdir) { this->localdir = localdir; }
	void setqsubdir(std::string qsubdir) { this->qsubdir = qsubdir; }
};


UserVariables::UserVariables(std::vector <std::string> lines)
{
	UserVariables::setmakenewfile(lines[2]);
	UserVariables::setfile(lines[4]);
	UserVariables::setlogin(lines[6]);
	UserVariables::sethpcdir(lines[8]);
	UserVariables::setlocaldir(lines[10]);
	UserVariables::setqsubdir(lines[12]);
};

class GeometryVariables : public Variables
{

private:
	int maxnumele;
	int numele;
	int maxstatevar;
	int statevar;
	int maxnodesperfe;
	int nodesperfe;

public:
	GeometryVariables(std::vector <std::string> lines);   // This is the constructor declaration
	~GeometryVariables() {};  // This is the destructor: declaration

	int getmaxnumele() { return this->maxnumele; }
	int getnumele() { return this->numele; }
	int getmaxstatevar() { return this->maxstatevar; }
	int getstatevar() { return this->statevar; }
	int getmaxnodesperfe() { return this->maxnodesperfe; }
	int getnodesperfe() { return this->nodesperfe; }
	
	void setmaxnumele(int maxnumele) { this->maxnumele = maxnumele; }
	void setnumele(int numele) { this->numele = numele; }
	void setmaxstatevar(int maxstatevar) { this->maxstatevar = maxstatevar; }
	void setstatevar(int statevar) { this->statevar = statevar; }
	void setmaxnodesperfe(int maxnodesperfe) { this->maxnodesperfe = maxnodesperfe; }
	void setnodesperfe(int nodesperfe) { this->nodesperfe = nodesperfe; }
};

GeometryVariables::GeometryVariables(std::vector <std::string> lines)
{
	GeometryVariables::setmaxnumele(1930);
	GeometryVariables::setnumele(190);
	GeometryVariables::setmaxstatevar(2);
	GeometryVariables::setstatevar(2);
	
};

class Inputfile
{
public:
	std::vector <std::string> readgeometry();
	std::vector <std::string> readuser();
	std::vector <std::string> readproblem();
	std::vector <std::string> readbonus();
	std::vector <std::string> readmaterial(std::string material);
};

std::vector <std::string> Inputfile::readgeometry()
{
	std::string line;
	std::vector <std::string> lines;
	std::ifstream infile;

	infile.open("geometry.txt");
	if (infile.is_open())
	{
		while (getline(infile, line)) {
			std::cout << line << std::endl;
			lines.push_back(line);
		}
		infile.close();
	}
	else
	{
		std::cout << "Unable to open user.txt" << std::endl;
	}

	return lines;
}

std::vector <std::string> Inputfile::readuser() 
{
	std::string line;
	std::vector <std::string> lines;
	std::ifstream infile;

	infile.open("user.txt");
	if (infile.is_open()) 
	{
		while (getline(infile, line)) {
			std::cout << line << std::endl;
			lines.push_back(line);
		}
		infile.close();
	}
	else 
	{
		std::cout << "Unable to open user.txt" << std::endl;
	}
	
	return lines;
}

std::vector <std::string> Inputfile::readproblem()
{
	std::string line;
	std::vector <std::string> lines;
	std::ifstream infile;

	infile.open("problemdata.txt");
	if (infile.is_open())
	{
		while (getline(infile, line)) {
			std::cout << line << std::endl;
			lines.push_back(line);
		}
		infile.close();
	}
	else
	{
		std::cout << "Unable to open user.txt" << std::endl;
	}

	return lines;
};

std::vector <std::string> Inputfile::readbonus()
{
	std::string line;
	std::vector <std::string> lines;
	std::ifstream infile;

	infile.open("bonusvariables.txt");
	if (infile.is_open())
	{
		while (getline(infile, line)) {
			std::cout << line << std::endl;
			lines.push_back(line);
		}
		infile.close();
	}
	else
	{
		std::cout << "Unable to open user.txt" << std::endl;
	}

	return lines;
};

std::vector <std::string> Inputfile::readmaterial(std::string material)
{
	std::string line;
	std::vector <std::string> lines;
	std::ifstream infile;

	infile.open(material +".txt");
	if (infile.is_open())
	{
		while (getline(infile, line)) {
			std::cout << line << std::endl;
			lines.push_back(line);
		}
		infile.close();
	}
	else
	{
		std::cout << "Unable to open user.txt" << std::endl;
	}

	return lines;
};

class Outfile
{

public:
	std::ofstream outfile;
	void writecontrol(UserVariables User, ProblemVariables Problem, BonusVariables Bonus, Material Steel, Material Pvb, Material Glas, Material Rock);
	void writeelements(GeometryVariables Geometry, UserVariables User, ProblemVariables Problem, BonusVariables Bonus, Material Steel, Material Pvb, Material Glas, Material Rock);
	void writeinteractions(UserVariables User, ProblemVariables Problem, BonusVariables Bonus, Material Steel, Material Pvb, Material Glas, Material Rock);
	void writenodes(UserVariables User, ProblemVariables Problem, BonusVariables Bonus, Material Steel, Material Pvb, Material Glas, Material Rock);
	void writeproperties(UserVariables User, ProblemVariables Problem, BonusVariables Bonus, Material Steel, Material Pvb, Material Glas, Material Rock);
	void writeoutput(UserVariables User, ProblemVariables Problem, BonusVariables Bonus, Material Steel, Material Pvb, Material Glas, Material Rock);
};

void Outfile::writecontrol(UserVariables User, ProblemVariables Problem, BonusVariables Bonus, Material Steel, Material Pvb, Material Glas, Material Rock)
{
	std::ofstream outfile;
	outfile.open(User.getfile(), std::ios::out | std::ios::trunc | std::ios::binary); ;
	if (outfile.is_open())
	{
		outfile << "\t/*   Control     */\n";
		outfile << "\t/YD/YDC/MCSTEP "		 												<< Problem.getmaxndt()					<< "\n"; //Maximum number of steps
		outfile << "\t\t/YD/YDC/NCSTEP "													<< Problem.getcurrentdt()				<< "\n"; //Current / Actual number of time steps (cannot be greater than /YD/YDC/MCSTEP)
		outfile << "\t\t/YD/YDC/DCGRAY "	<< std::setprecision(11)	<< std::scientific	<< Problem.getgravity()					<< "\n"; //Acceleration of gravity (in y direction)
		outfile << "\t\t/YD/YDC/DCSIZC "	<< std::setprecision(11)	<< std::scientific	<< Problem.getcoordsize()				<< "\n"; //Maximum size of coordinate in any direction (size of physical space – corresponding outputs are normalized using this value)
		outfile << "\t\t/YD/YDC/DCSIZF "	<< std::setprecision(11)	<< std::scientific	<< Problem.getmaxforce()				<< "\n"; //Maximum size of force in any direction(corresponding outputs are normalized using this value)
		outfile << "\t\t/YD/YDC/DCSIZS "	<< std::setprecision(11)	<< std::scientific	<< Problem.getmaxstress()				<< "\n"; //Maximum stress (buffer size)
		outfile << "\t\t/YD/YDC/DCSIZV "	<< std::setprecision(11)	<< std::scientific	<< Problem.getmaxv()					<< "\n"; //Maximum size of velocity in any direction (corresponding outputs are normalized using this value)
		outfile << "\t\t/YD/YDC/DCSIZD "	<< std::setprecision(11)	<< std::scientific	<< Problem.getmaxdispl()				<< "\n"; //Maximum size of displacement in any direction
		outfile << "\t\t/YD/YDC/DCSIZA "	<< std::setprecision(11)	<< std::scientific	<< Problem.getmaxjointaperture()		<< "\n"; //Maximum joint aperture
		outfile << "\t\t/YD/YDC/DCSTEC "	<< std::setprecision(11)	<< std::scientific	<< Problem.getdt()						<< "\n"; //Size of the time step (see in the book how to calculate critical-maximum time step)
		outfile << "\t\t/YD/YDC/DCTIME "	<< std::setprecision(11)	<< std::scientific	<< Problem.getcurrentt()				<< "\n"; //Current time. i.e. time at start of this run.
		outfile << "\t\t/YD/YDC/DCRMPT "	<< std::setprecision(11)	<< std::scientific	<< Problem.getloadrampingstage()		<< "\n"; //Load ramping stage
		outfile << "\t\t/YD/YDC/DCGRST "	<< std::setprecision(11)	<< std::scientific	<< Problem.getgravitysettlingstage()	<< "\n"; //Gravity settling stage
		outfile << "\t\t/YD/YDC/ICOUTF "													<< Problem.getoutfreq()					<< "\n"; //Output frequency – every so many time steps complete state of the system is recorded in a file with extension .ym which can be visualized using M program, which is FEM/DEM Visualizer accompanying Y program.
		outfile << "\t\t/YD/YDC/ICOUTI "													<< 0									<< "\n"; //Current number of iterations; Serial number of first output associated with this run.
		outfile << "\t\t/YD/YDC/ICSAVF "													<< Problem.getrestartfreq()				<< "\n"; //Restart save frequency
		outfile << "\t\t/YD/YDC/ICOUTP "													<< 4									<< "\n"; //Number of characters for each number in output file (for example: three characters is equivalent to six significant digits)
		outfile << "\t\t/YD/YDC/ICFMTY "													<< 0									<< "\n"; //
		outfile << "\t\t/YD/YDC/ICIATY "													<< 0									<< "\n"; //
		outfile.close();
	}
	else
	{
		std::cout << "Unable to open " << User.getfile() << std::endl;
	};
};

void Outfile::writeelements(GeometryVariables Geometry, UserVariables User, ProblemVariables Problem, BonusVariables Bonus, Material Steel, Material Pvb, Material Glas, Material Rock)
{
	std::ofstream outfile;
	outfile.open(User.getfile(), std::ios::out | std::ios::app | std::ios::binary); ;
	if (outfile.is_open())
	{
		outfile << "\t/*   Elements     */\n"; 
		outfile << "\t/YD/YDE/MELEM " << Geometry.getmaxnumele() << "\n";						//Maximum number of finite elements(with fracture the actual number of finite elements increases during the run, however, it should not exceed this number)
		outfile << "\t/YD/YDE/NELEM " << Geometry.getnumele() << "\n";							//Actual number of finite elements at the beginning of this run.
		outfile << "\t/YD/YDE/MELST " << Geometry.getmaxstatevar() << "\n";						//Maximum number of state variables per finite element.
		outfile << "\t/YD/YDE/NELST " << Geometry.getstatevar() << "\n";						//Actual number of state variables per finite element.
		outfile << "\t/YD/YDE/MELNO " << Geometry.getmaxnodesperfe() << "\n";					//Maximum number of nodes per finite element.
		outfile << "\t/YD/YDE/NELNO " << Geometry.getnodesperfe() << "\n";						//Actual number of nodes per finite element.
		outfile << "\t/YD/YDE/D2ELST " << 21 << " " << Geometry.getnumele() << " " << 0 << "\n";//[MELST][MELEM] array containing state variables for all finite elements.
		outfile << "\t/YD/YDE/I1ELCF " << Geometry.getnumele() << "\n";							//Head of a list of contacting couples for each finite element.
		outfile << "\t";
		for (int i = 0; i < Geometry.getnumele(); i++) { outfile << "-1 "; };
		outfile << "\n";
		outfile << "\t/YD/YDE/I1ELTY " << Geometry.getnumele() << "\n";
		outfile << "\t";
		for (int i = 0; i < Geometry.getnumele(); i++) { outfile << "-1 "; };
		outfile << "\n";
		outfile << "\t/YD/YDE/I1ELPR " << Geometry.getnumele() << "\n"; //Set of properties associated with each element.
		outfile << "\t1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0 	\n";
		outfile << "\t/YD/YDE/I2ELTO    21     3   190\n";//[MELNO][MELEM] topology array containing nodes for each finite element.
		outfile << "\t4      9      7\n";
		outfile << "\t130    129    12\n";
		outfile << "\t0      1      2\n";
		outfile << "\t128    127    12\n";
		outfile << "\t127    126    12\n";
		outfile << "\t129    128    12\n";
		outfile << "\t9     14     16\n";
		outfile << "\t1      4      3\n";
		outfile << "\t7      9     16\n";
		outfile << "\t124    128    12\n";
		outfile << "\t124    123    11\n";
		outfile << "\t118    123    11\n";
		outfile << "\t124    118    12\n";
		outfile << "\t123    127    12\n";
		outfile << "\t3      4      7\n";
		outfile << "\t3      7      8\n";
		outfile << "\t3      8      5\n";
		outfile << "\t8      7     12\n";
		outfile << "\t129    124    12\n";
		outfile << "\t1      3      2\n";
		outfile << "\t12      7     16\n";
		outfile << "\t5      8     10\n";
		outfile << "\t121    118    11\n";
		outfile << "\t119    123    12\n";
		outfile << "\t3      5      2\n";
		outfile << "\t8     12     13\n";
		outfile << "\t124    121    12\n";
		outfile << "\t118    119    11\n";
		outfile << "\t13     18     22\n";
		outfile << "\t115    117    11\n";
		outfile << "\t116    118    11\n";
		outfile << "\t10      8     13\n";
		outfile << "\t121    116    12\n";
		outfile << "\t5     10      6\n";
		outfile << "\t17     11     15\n";
		outfile << "\t15     11     10\n";
		outfile << "\t17     15     21\n";
		outfile << "\t17     21     23\n";
		outfile << "\t21     15     19\n";
		outfile << "\t19     15     10\n";
		outfile << "\t21     19     25\n";
		outfile << "\t25     19     24\n";
		outfile << "\t24     19     20\n";
		outfile << "\t24     20     26\n";
		outfile << "\t20     19     13\n";
		outfile << "\t25     24     29\n";
		outfile << "\t29     24     30\n";
		outfile << "\t25     29     31\n";
		outfile << "\t25     31     27\n";
		outfile << "\t25     27     21\n";
		outfile << "\t21     27     23\n";
		outfile << "\t23     27     28\n";
		outfile << "\t27     31     32\n";
		outfile << "\t32     31     36\n";
		outfile << "\t36     31     33\n";
		outfile << "\t36     33     39\n";
		outfile << "\t39     33     38\n";
		outfile << "\t32     36     37\n";
		outfile << "\t37     36     41\n";
		outfile << "\t37     41     42\n";
		outfile << "\t42     41     50\n";
		outfile << "\t42     50     55\n";
		outfile << "\t55     50     66\n";
		outfile << "\t66     50     62\n";
		outfile << "\t62     50     46\n";
		outfile << "\t46     50     41\n";
		outfile << "\t46     41     39\n";
		outfile << "\t46     39     45\n";
		outfile << "\t39     41     36\n";
		outfile << "\t62     46     58\n";
		outfile << "\t31     29     33\n";
		outfile << "\t33     29     34\n";
		outfile << "\t55     66     70\n";
		outfile << "\t55     70     64\n";
		outfile << "\t64     70     77\n";
		outfile << "\t70     66     78\n";
		outfile << "\t70     78     83\n";
		outfile << "\t83     78     94\n";
		outfile << "\t83     94     99\n";
		outfile << "\t99     94    103\n";
		outfile << "\t103     94    10\n";
		outfile << "\t101     94     9\n";
		outfile << "\t101     90    10\n";
		outfile << "\t90     94     78\n";
		outfile << "\t90     78     75\n";
		outfile << "\t75     78     66\n";
		outfile << "\t75     66     62\n";
		outfile << "\t75     62     72\n";
		outfile << "\t90     75     87\n";
		outfile << "\t99    103    104\n";
		outfile << "\t104    103    10\n";
		outfile << "\t104    108    10\n";
		outfile << "\t109    108    11\n";
		outfile << "\t99    104    102\n";
		outfile << "\t102    104    10\n";
		outfile << "\t108    103    10\n";
		outfile << "\t106    103    10\n";
		outfile << "\t106    101    10\n";
		outfile << "\t108    106    11\n";
		outfile << "\t111    106    11\n";
		outfile << "\t108    111    11\n";
		outfile << "\t112    111    11\n";
		outfile << "\t112    116    11\n";
		outfile << "\t112    115    10\n";
		outfile << "\t104    109    10\n";
		outfile << "\t107    109    11\n";
		outfile << "\t70     83     77\n";
		outfile << "\t77     83     93\n";
		outfile << "\t32     37     35\n";
		outfile << "\t35     37     40\n";
		outfile << "\t83     99     93\n";
		outfile << "\t93     99    102\n";
		outfile << "\t42     55     49\n";
		outfile << "\t49     55     64\n";
		outfile << "\t27     32     28\n";
		outfile << "\t28     32     35\n";
		outfile << "\t37     42     40\n";
		outfile << "\t40     42     49\n";
		outfile << "\t120    125    12\n";
		outfile << "\t120    116    11\n";
		outfile << "\t45     58     46\n";
		outfile << "\t26     30     24\n";
		outfile << "\t100    105    10\n";
		outfile << "\t117    115    11\n";
		outfile << "\t117    119    12\n";
		outfile << "\t18     13     12\n";
		outfile << "\t18     12     16\n";
		outfile << "\t6      2      5 \n";
		outfile << "\t6     10     11 \n";
		outfile << "\t72     87     75\n";
		outfile << "\t34     38     33\n";
		outfile << "\t110    114    11\n";
		outfile << "\t58     72     62\n";
		outfile << "\t105    110    10\n";
		outfile << "\t30     34     29\n";
		outfile << "\t22     26     20\n";
		outfile << "\t87    100     90\n";
		outfile << "\t38     45     39\n";
		outfile << "\t10     13     19\n";
		outfile << "\t20     13     22\n";
		outfile << "\t113    109    11\n";
		outfile << "\t116    111    11\n";
		outfile << "\t96     92     89\n";
		outfile << "\t59     52     61\n";
		outfile << "\t65     71     68\n";
		outfile << "\t43     47     53\n";
		outfile << "\t73     67     69\n";
		outfile << "\t91     95     85\n";
		outfile << "\t79     84     80\n";
		outfile << "\t88     82     81\n";
		outfile << "\t48     44     54\n";
		outfile << "\t51     56     57\n";
		outfile << "\t97     98     86\n";
		outfile << "\t52     48     61\n";
		outfile << "\t92     88     81\n";
		outfile << "\t71     79     80\n";
		outfile << "\t47     51     53\n";
		outfile << "\t95     97     85\n";
		outfile << "\t67     59     61\n";
		outfile << "\t84     91     80\n";
		outfile << "\t44     43     53\n";
		outfile << "\t82     73     81\n";
		outfile << "\t98     96     89\n";
		outfile << "\t56     65     68\n";
		outfile << "\t89     92     81\n";
		outfile << "\t54     44     53\n";
		outfile << "\t69     67     61\n";
		outfile << "\t68     71     80\n";
		outfile << "\t86     98     89\n";
		outfile << "\t57     56     68\n";
		outfile << "\t48     54     61\n";
		outfile << "\t73     69     81\n";
		outfile << "\t97     86     85\n";
		outfile << "\t91     85     80\n";
		outfile << "\t51     57     53\n";
		outfile << "\t89     81     76\n";
		outfile << "\t86     89     76\n";
		outfile << "\t85     86     74\n";
		outfile << "\t61     54     63\n";
		outfile << "\t53     57     63\n";
		outfile << "\t69     61     63\n";
		outfile << "\t80     85     74\n";
		outfile << "\t54     53     63\n";
		outfile << "\t57     68     74\n";
		outfile << "\t81     69     76\n";
		outfile << "\t68     80     74\n";
		outfile << "\t74     86     76\n";
		outfile << "\t74     76     63\n";
		outfile << "\t63     76     69\n";
		outfile << "\t74     63     57\n";
		outfile << "\t/YD/YDE/I2ELJP   21  3  0\n"; 
		outfile.close();
	}
	else
	{
		std::cout << "Unable to open " << User.getfile() << std::endl;
	}
}

void Outfile::writeinteractions(UserVariables User, ProblemVariables Problem, BonusVariables Bonus, Material Steel, Material Pvb, Material Glas, Material Rock)
{
	std::ofstream outfile;
	outfile.open(User.getfile(), std::ios::out | std::ios::app | std::ios::binary); ;
	if (outfile.is_open())
	{
		outfile << "\t/*   Interactions     */\n";
		outfile << "\t/YD/YDI/MICOUP " << Problem.getmaxcontcouples() << "\n"; //Maximum number of contacting couples of finite elements.
		outfile << "\t/YD/YDI/NICOUP " << Problem.getcontcouples()	   << "\n"; //Actual number of contacting couples of finite elements (always set to zero)
		outfile << "\t/YD/YDI/IIECFF -2\n";//Internal variable used for contact (always set to -2)
		outfile << "\t/YD/YDI/DIEDI +2e+002\n"; //Internal variable (travel since last detection) used to trigger contact detection (always set to a large number in order to trigger contact detection immediately at the start of this run)
		outfile << "\t/YD/YDI/DIEZON " << std::scientific << Problem.getbuffersize() << "\n"; //Contact detection buffer size, Size of the buffer around each finite element for contact detection purposes (usually 1/5 of the size of the smallest finite element)
		outfile << "\t/YD/YDI/D1IESL 0\n"; //[MELEM] array used to store sliding distance between couples in contact
		outfile << "\t/YD/YDI/I1IECN 0\n"; //[MICOUP] array used to store next contacting couple in the list.
		outfile << "\t/YD/YDI/I1IECT 0\n"; //[MICOUP] array used to store target finite element for each contacting couple.
	}
	else
	{
		std::cout << "Unable to open " << User.getfile() << std::endl;
	}
}

void Outfile::writenodes(UserVariables User, ProblemVariables Problem, BonusVariables Bonus, Material Steel, Material Pvb, Material Glas, Material Rock)
{
	std::ofstream outfile;
	outfile.open(User.getfile(), std::ios::out | std::ios::app | std::ios::binary); ;
	if (outfile.is_open())
	{
	outfile << "/*   Nodes     */\n";
	outfile << "/YD/YDN/MNODIM 3\n";	//Maximum number of degrees of freedom per node (usually 2 in 2D)
	outfile << "/YD/YDN/NNODIM 2\n";	//Actual number of degrees of freedom per node (usually 2 in 2D)
	outfile << "/YD/YDN/MNOPO 6600\n";	//Maximum number of nodes.
	outfile << "/YD/YDN/NNOPO 131\n";	//Actual number of nodes.
	outfile << "/YD/YDN/D2NCC 21 2 131\n"; //[MNODIM][MNOPO] array containing the current coordinates of the nodes.
	outfile << "+2.000000000000e+000 +1.000000000000e+001\n";
	outfile << "+1.500000000000e+000 +1.000000000000e+001\n";
	outfile << "+2.000000000000e+000 +9.444444444400e+000\n";
	outfile << "+1.250000000000e+000 +9.523686027900e+000\n";
	outfile << "+1.000000000000e+000 +1.000000000000e+001\n";
	outfile << "+1.507600006800e+000 +9.096814324300e+000\n";
	outfile << "+2.000000000000e+000 +8.888888888900e+000\n";
	outfile << "+7.500000000000e-001 +9.523686027900e+000\n";
	outfile << "+1.000000000000e+000 +9.047372055800e+000\n";
	outfile << "+5.000000000000e-001 +1.000000000000e+001\n";
	outfile << "+1.288000033900e+000 +8.579680204300e+000\n";
	outfile << "+2.000000000000e+000 +8.333333333300e+000\n";
	outfile << "+4.923999932200e-001 +9.096814324300e+000\n";
	outfile << "+7.119999660800e-001 +8.579680204300e+000\n";
	outfile << "+0.000000000000e+000 +1.000000000000e+001\n";
	outfile << "+1.523686027900e+000 +8.055555555600e+000\n";
	outfile << "+0.000000000000e+000 +9.444444444400e+000\n";
	outfile << "+2.000000000000e+000 +7.777777777800e+000\n";
	outfile << "+0.000000000000e+000 +8.888888888900e+000\n";
	outfile << "+1.017841164100e+000 +7.922412906400e+000\n";
	outfile << "+4.655371738300e-001 +7.996262427700e+000\n";
	outfile << "+1.521888636600e+000 +7.505209625400e+000\n";
	outfile << "+0.000000000000e+000 +8.333333333300e+000\n";
	outfile << "+2.000000000000e+000 +7.222222222200e+000\n";
	outfile << "+5.692100026200e-001 +7.508582632800e+000\n";
	outfile << "+1.044566307600e+000 +7.231919694600e+000\n";
	outfile << "+0.000000000000e+000 +7.777777777800e+000\n";
	outfile << "+1.519916414700e+000 +6.955256778100e+000\n";
	outfile << "+2.000000000000e+000 +6.666666666700e+000\n";
	outfile << "+5.672933977200e-001 +6.958584609000e+000\n";
	outfile << "+0.000000000000e+000 +7.222222222200e+000\n";
	outfile << "+1.042644393800e+000 +6.681922711800e+000\n";
	outfile << "+1.517994612700e+000 +6.405260114400e+000\n";
	outfile << "+5.653724697700e-001 +6.408588645100e+000\n";
	outfile << "+0.000000000000e+000 +6.666666666700e+000\n";
	outfile << "+2.000000000000e+000 +6.111111111100e+000\n";
	outfile << "+1.040722674900e+000 +6.131926063700e+000\n";
	outfile << "+1.516072873100e+000 +5.855263470400e+000\n";
	outfile << "+0.000000000000e+000 +6.111111111100e+000\n";
	outfile << "+5.634507491500e-001 +5.858592002100e+000\n";
	outfile << "+2.000000000000e+000 +5.555555555600e+000\n";
	outfile << "+1.038800940400e+000 +5.581929420800e+000\n";
	outfile << "+1.514151137300e+000 +5.305266827700e+000\n";
	outfile << "-5.845849858500e-001 +5.909631994800e+000\n";
	outfile << "-3.451392653500e-001 +5.755749573700e+000\n";
	outfile << "+0.000000000000e+000 +5.555555555600e+000\n";
	outfile << "+5.615290082200e-001 +5.308595371200e+000\n";
	outfile << "-8.576851601600e-001 +5.989821441700e+000\n";
	outfile << "-1.587464668300e-001 +5.540640816900e+000\n";
	outfile << "+2.000000000000e+000 +5.000000000000e+000\n";
	outfile << "+1.036879205000e+000 +5.031932778100e+000\n";
	outfile << "-1.142314836400e+000 +5.989821442200e+000\n";
	outfile << "-4.050702629600e-002 +5.281732556500e+000\n";
	outfile << "-8.044148338100e-001 +5.666101901100e+000\n";
	outfile << "-4.753414119500e-001 +5.454619253700e+000\n";
	outfile << "+1.512229401800e+000 +4.755270185000e+000\n";
	outfile << "-1.415415011800e+000 +5.909631995900e+000\n";
	outfile << "-1.195585163900e+000 +5.666101900700e+000\n";
	outfile << "+0.000000000000e+000 +5.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +5.000000000000e+000\n";
	outfile << "-4.943403017200e-002 +5.000000000000e+000\n";
	outfile << "-3.128433843500e-001 +5.098798205800e+000\n";
	outfile << "+5.596072727300e-001 +4.758598728500e+000\n";
	outfile << "-7.675007363000e-001 +5.188298771100e+000\n";
	outfile << "+2.000000000000e+000 +4.444444444400e+000\n";
	outfile << "-1.654860733800e+000 +5.755749574500e+000\n";
	outfile << "+1.034957469500e+000 +4.481936135400e+000\n";
	outfile << "-4.050702693800e-002 +4.718267441300e+000\n";
	outfile << "-1.524658587500e+000 +5.454619254300e+000\n";
	outfile << "-3.685127272600e-001 +4.711609424200e+000\n";
	outfile << "+1.510307666300e+000 +4.205273542300e+000\n";
	outfile << "-1.841253532500e+000 +5.540640817900e+000\n";
	outfile << "+0.000000000000e+000 +4.444444444400e+000\n";
	outfile << "-1.587464680800e-001 +4.459359181100e+000\n";
	outfile << "-1.301051933900e+000 +5.043729668000e+000\n";
	outfile << "+5.576855372300e-001 +4.208602085900e+000\n";
	outfile << "-9.147556989700e-001 +4.677131044300e+000\n";
	outfile << "+2.000000000000e+000 +3.888888888900e+000\n";
	outfile << "+1.033035734000e+000 +3.931939492800e+000\n";
	outfile << "-1.959492973400e+000 +5.281732557400e+000\n";
	outfile << "-1.745348182200e+000 +5.107228452900e+000\n";
	outfile << "-6.246748171100e-001 +4.415982614100e+000\n";
	outfile << "-3.451392671000e-001 +4.244250424700e+000\n";
	outfile << "+1.508385930800e+000 +3.655276899700e+000\n";
	outfile << "-2.000000000000e+000 +5.000000001200e+000\n";
	outfile << "-1.631487273600e+000 +4.711609426600e+000\n";
	outfile << "-1.329248315600e+000 +4.508721584400e+000\n";
	outfile << "+0.000000000000e+000 +3.888888888900e+000\n";
	outfile << "-5.845849879600e-001 +4.090368004200e+000\n";
	outfile << "-1.057919239000e+000 +4.285240936100e+000\n";
	outfile << "+5.557638017400e-001 +3.658605443200e+000\n";
	outfile << "-1.959492974100e+000 +4.718267444800e+000\n";
	outfile << "-8.576851624600e-001 +4.010178558000e+000\n";
	outfile << "+2.000000000000e+000 +3.333333333300e+000\n";
	outfile << "+1.031113998500e+000 +3.381942850100e+000\n";
	outfile << "-1.841253533800e+000 +4.459359184100e+000\n";
	outfile << "-1.142314838700e+000 +4.010178558200e+000\n";
	outfile << "-1.654860735300e+000 +4.244250426900e+000\n";
	outfile << "-1.415415013100e+000 +4.090368004700e+000\n";
	outfile << "+1.506464195300e+000 +3.105280257000e+000\n";
	outfile << "+0.000000000000e+000 +3.333333333300e+000\n";
	outfile << "+5.538420662400e-001 +3.108608800500e+000\n";
	outfile << "+2.000000000000e+000 +2.777777777800e+000\n";
	outfile << "+1.029192263000e+000 +2.831946207500e+000\n";
	outfile << "+1.504542459800e+000 +2.555283614400e+000\n";
	outfile << "+0.000000000000e+000 +2.777777777800e+000\n";
	outfile << "+5.519203307400e-001 +2.558612157900e+000\n";
	outfile << "+2.000000000000e+000 +2.222222222200e+000\n";
	outfile << "+1.027270527500e+000 +2.281949564800e+000\n";
	outfile << "+1.502620724300e+000 +2.005286971700e+000\n";
	outfile << "+0.000000000000e+000 +2.222222222200e+000\n";
	outfile << "+5.499985952500e-001 +2.008615515200e+000\n";
	outfile << "+1.025348792000e+000 +1.731952922100e+000\n";
	outfile << "+2.000000000000e+000 +1.666666666700e+000\n";
	outfile << "+0.000000000000e+000 +1.666666666700e+000\n";
	outfile << "+1.288000033900e+000 +1.420319795700e+000\n";
	outfile << "+7.119999660800e-001 +1.420319795700e+000\n";
	outfile << "+2.000000000000e+000 +1.111111111100e+000\n";
	outfile << "+1.000000000000e+000 +9.526279441600e-001\n";
	outfile << "+1.507600006800e+000 +9.031856757100e-001\n";
	outfile << "+0.000000000000e+000 +1.111111111100e+000\n";
	outfile << "+4.923999932200e-001 +9.031856757100e-001\n";
	outfile << "+2.000000000000e+000 +5.555555555600e-001\n";
	outfile << "+1.250000000000e+000 +4.763139720800e-001\n";
	outfile << "+7.500000000000e-001 +4.763139720800e-001\n";
	outfile << "+0.000000000000e+000 +5.555555555600e-001\n";
	outfile << "+2.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+1.500000000000e+000 +0.000000000000e+000\n";
	outfile << "+1.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e-001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n\n";
	outfile << "/YD/YDN/D2NCI 21 2 131\n"; //[MNODIM][MNOPO] array containing the initial coordinates of the nodes(usually corresponding to undeformed system)
	outfile << "+2.000000000000e+000 +1.000000000000e+001\n";
	outfile << "+1.500000000000e+000 +1.000000000000e+001\n";
	outfile << "+2.000000000000e+000 +9.444444444400e+000\n";
	outfile << "+1.250000000000e+000 +9.523686027900e+000\n";
	outfile << "+1.000000000000e+000 +1.000000000000e+001\n";
	outfile << "+1.507600006800e+000 +9.096814324300e+000\n";
	outfile << "+2.000000000000e+000 +8.888888888900e+000\n";
	outfile << "+7.500000000000e-001 +9.523686027900e+000\n";
	outfile << "+1.000000000000e+000 +9.047372055800e+000\n";
	outfile << "+5.000000000000e-001 +1.000000000000e+001\n";
	outfile << "+1.288000033900e+000 +8.579680204300e+000\n";
	outfile << "+2.000000000000e+000 +8.333333333300e+000\n";
	outfile << "+4.923999932200e-001 +9.096814324300e+000\n";
	outfile << "+7.119999660800e-001 +8.579680204300e+000\n";
	outfile << "+0.000000000000e+000 +1.000000000000e+001\n";
	outfile << "+1.523686027900e+000 +8.055555555600e+000\n";
	outfile << "+0.000000000000e+000 +9.444444444400e+000\n";
	outfile << "+2.000000000000e+000 +7.777777777800e+000\n";
	outfile << "+0.000000000000e+000 +8.888888888900e+000\n";
	outfile << "+1.017841164100e+000 +7.922412906400e+000\n";
	outfile << "+4.655371738300e-001 +7.996262427700e+000\n";
	outfile << "+1.521888636600e+000 +7.505209625400e+000\n";
	outfile << "+0.000000000000e+000 +8.333333333300e+000\n";
	outfile << "+2.000000000000e+000 +7.222222222200e+000\n";
	outfile << "+5.692100026200e-001 +7.508582632800e+000\n";
	outfile << "+1.044566307600e+000 +7.231919694600e+000\n";
	outfile << "+0.000000000000e+000 +7.777777777800e+000\n";
	outfile << "+1.519916414700e+000 +6.955256778100e+000\n";
	outfile << "+2.000000000000e+000 +6.666666666700e+000\n";
	outfile << "+5.672933977200e-001 +6.958584609000e+000\n";
	outfile << "+0.000000000000e+000 +7.222222222200e+000\n";
	outfile << "+1.042644393800e+000 +6.681922711800e+000\n";
	outfile << "+1.517994612700e+000 +6.405260114400e+000\n";
	outfile << "+5.653724697700e-001 +6.408588645100e+000\n";
	outfile << "+0.000000000000e+000 +6.666666666700e+000\n";
	outfile << "+2.000000000000e+000 +6.111111111100e+000\n";
	outfile << "+1.040722674900e+000 +6.131926063700e+000\n";
	outfile << "+1.516072873100e+000 +5.855263470400e+000\n";
	outfile << "+0.000000000000e+000 +6.111111111100e+000\n";
	outfile << "+5.634507491500e-001 +5.858592002100e+000\n";
	outfile << "+2.000000000000e+000 +5.555555555600e+000\n";
	outfile << "+1.038800940400e+000 +5.581929420800e+000\n";
	outfile << "+1.514151137300e+000 +5.305266827700e+000\n";
	outfile << "-5.845849858500e-001 +5.909631994800e+000\n";
	outfile << "-3.451392653500e-001 +5.755749573700e+000\n";
	outfile << "+0.000000000000e+000 +5.555555555600e+000\n";
	outfile << "+5.615290082200e-001 +5.308595371200e+000\n";
	outfile << "-8.576851601600e-001 +5.989821441700e+000\n";
	outfile << "-1.587464668300e-001 +5.540640816900e+000\n";
	outfile << "+2.000000000000e+000 +5.000000000000e+000\n";
	outfile << "+1.036879205000e+000 +5.031932778100e+000\n";
	outfile << "-1.142314836400e+000 +5.989821442200e+000\n";
	outfile << "-4.050702629600e-002 +5.281732556500e+000\n";
	outfile << "-8.044148338100e-001 +5.666101901100e+000\n";
	outfile << "-4.753414119500e-001 +5.454619253700e+000\n";
	outfile << "+1.512229401800e+000 +4.755270185000e+000\n";
	outfile << "-1.415415011800e+000 +5.909631995900e+000\n";
	outfile << "-1.195585163900e+000 +5.666101900700e+000\n";
	outfile << "+0.000000000000e+000 +5.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +5.000000000000e+000\n";
	outfile << "-4.943403017200e-002 +5.000000000000e+000\n";
	outfile << "-3.128433843500e-001 +5.098798205800e+000\n";
	outfile << "+5.596072727300e-001 +4.758598728500e+000\n";
	outfile << "-7.675007363000e-001 +5.188298771100e+000\n";
	outfile << "+2.000000000000e+000 +4.444444444400e+000\n";
	outfile << "-1.654860733800e+000 +5.755749574500e+000\n";
	outfile << "+1.034957469500e+000 +4.481936135400e+000\n";
	outfile << "-4.050702693800e-002 +4.718267441300e+000\n";
	outfile << "-1.524658587500e+000 +5.454619254300e+000\n";
	outfile << "-3.685127272600e-001 +4.711609424200e+000\n";
	outfile << "+1.510307666300e+000 +4.205273542300e+000\n";
	outfile << "-1.841253532500e+000 +5.540640817900e+000\n";
	outfile << "+0.000000000000e+000 +4.444444444400e+000\n";
	outfile << "-1.587464680800e-001 +4.459359181100e+000\n";
	outfile << "-1.301051933900e+000 +5.043729668000e+000\n";
	outfile << "+5.576855372300e-001 +4.208602085900e+000\n";
	outfile << "-9.147556989700e-001 +4.677131044300e+000\n";
	outfile << "+2.000000000000e+000 +3.888888888900e+000\n";
	outfile << "+1.033035734000e+000 +3.931939492800e+000\n";
	outfile << "-1.959492973400e+000 +5.281732557400e+000\n";
	outfile << "-1.745348182200e+000 +5.107228452900e+000\n";
	outfile << "-6.246748171100e-001 +4.415982614100e+000\n";
	outfile << "-3.451392671000e-001 +4.244250424700e+000\n";
	outfile << "+1.508385930800e+000 +3.655276899700e+000\n";
	outfile << "-2.000000000000e+000 +5.000000001200e+000\n";
	outfile << "-1.631487273600e+000 +4.711609426600e+000\n";
	outfile << "-1.329248315600e+000 +4.508721584400e+000\n";
	outfile << "+0.000000000000e+000 +3.888888888900e+000\n";
	outfile << "-5.845849879600e-001 +4.090368004200e+000\n";
	outfile << "-1.057919239000e+000 +4.285240936100e+000\n";
	outfile << "+5.557638017400e-001 +3.658605443200e+000\n";
	outfile << "-1.959492974100e+000 +4.718267444800e+000\n";
	outfile << "-8.576851624600e-001 +4.010178558000e+000\n";
	outfile << "+2.000000000000e+000 +3.333333333300e+000\n";
	outfile << "+1.031113998500e+000 +3.381942850100e+000\n";
	outfile << "-1.841253533800e+000 +4.459359184100e+000\n";
	outfile << "-1.142314838700e+000 +4.010178558200e+000\n";
	outfile << "-1.654860735300e+000 +4.244250426900e+000\n";
	outfile << "-1.415415013100e+000 +4.090368004700e+000\n";
	outfile << "+1.506464195300e+000 +3.105280257000e+000\n";
	outfile << "+0.000000000000e+000 +3.333333333300e+000\n";
	outfile << "+5.538420662400e-001 +3.108608800500e+000\n";
	outfile << "+2.000000000000e+000 +2.777777777800e+000\n";
	outfile << "+1.029192263000e+000 +2.831946207500e+000\n";
	outfile << "+1.504542459800e+000 +2.555283614400e+000\n";
	outfile << "+0.000000000000e+000 +2.777777777800e+000\n";
	outfile << "+5.519203307400e-001 +2.558612157900e+000\n";
	outfile << "+2.000000000000e+000 +2.222222222200e+000\n";
	outfile << "+1.027270527500e+000 +2.281949564800e+000\n";
	outfile << "+1.502620724300e+000 +2.005286971700e+000\n";
	outfile << "+0.000000000000e+000 +2.222222222200e+000\n";
	outfile << "+5.499985952500e-001 +2.008615515200e+000\n";
	outfile << "+1.025348792000e+000 +1.731952922100e+000\n";
	outfile << "+2.000000000000e+000 +1.666666666700e+000\n";
	outfile << "+0.000000000000e+000 +1.666666666700e+000\n";
	outfile << "+1.288000033900e+000 +1.420319795700e+000\n";
	outfile << "+7.119999660800e-001 +1.420319795700e+000\n";
	outfile << "+2.000000000000e+000 +1.111111111100e+000\n";
	outfile << "+1.000000000000e+000 +9.526279441600e-001\n";
	outfile << "+1.507600006800e+000 +9.031856757100e-001\n";
	outfile << "+0.000000000000e+000 +1.111111111100e+000\n";
	outfile << "+4.923999932200e-001 +9.031856757100e-001\n";
	outfile << "+2.000000000000e+000 +5.555555555600e-001\n";
	outfile << "+1.250000000000e+000 +4.763139720800e-001\n";
	outfile << "+7.500000000000e-001 +4.763139720800e-001\n";
	outfile << "+0.000000000000e+000 +5.555555555600e-001\n";
	outfile << "+2.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+1.500000000000e+000 +0.000000000000e+000\n";
	outfile << "+1.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e-001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n\n";
	outfile << "/YD/YDN/D2NFC 21 2 131\n"; //[MNODIM][MNOPO] array containing the current nodal forces.
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n\n";
	outfile << "/YD/YDN/D1NMCT 0\n\n";//[MNOPO] array containing the current mass for each node.
	outfile << "/YD/YDN/D2NVC 21 2 131\n";//[MNODIM][MNOPO] array containing the current velocities for each node.
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+5.000000000000e+001 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n";
	outfile << "+0.000000000000e+000 +0.000000000000e+000\n\n";
	outfile << "/YD/YDN/I1NOBF   131\n";//[MNOPO] array containing a flag indicating that a node is on the boundary (usually set to 1 for all nodes indicating the boundary. Any node with the flag set to zero is considered not to be on the boundary)
	outfile << "1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1	 \n";
	outfile << "/YD/YDN/I1NOPR   131\n";//[MNOPO] array containing an ID of a property set associated with each node.
	outfile << "4    4    3    3    4    3    3    3    3    4    3    3    3    3    4    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    3    4    4    4    4    4  \n";
	}
	else
	{
	std::cout << "Unable to open " << User.getfile() << std::endl;
	}
}

void Outfile::writeproperties(UserVariables User, ProblemVariables Problem, BonusVariables Bonus, Material Steel, Material Pvb, Material Glas, Material Rock)
{
	std::ofstream outfile;
	outfile.open(User.getfile(), std::ios::out | std::ios::app | std::ios::binary); ;
	if (outfile.is_open())
	{
		outfile << "/*   Properties     */\n";
		outfile << "/YD/YDP/MPROP 6\n";//Maximum number of properties sets.
		outfile << "/YD/YDP/NPROP 5\n";//Actual number of properties sets.
		outfile << "/YD/YDP/D1PCOH 5\n";//[MPROP] array containing internal cohesion
		outfile << "+8.00000000000e+006 +8.00000000000e+006 +8.00000000000e+006 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PICF 5\n";//[MPROP] array containing interface friction.
		outfile << "+6.00000000000e-001 +6.00000000000e-001 +6.00000000000e-001 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PEFR 5\n";//[MPROP] array containing sliding friction (from 0 – glass to 2.5 – rubber).
		outfile << "+6.00000000000e-001 +6.00000000000e-001 +6.00000000000e-001 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PEFT 5\n";//[MPROP] array containing tensile strength of material for finite elements for each of the properties.
		outfile << "+4.00000000000e+006 +4.00000000000e+006 +4.00000000000e+006 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PEGT 5\n"; //[MPROP] array containing mode 1 energy rate
		outfile << "+2.00000000000e+001 +2.00000000000e+001 +2.00000000000e+001 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PEGS 5\n"; //[MPROP] array containing damping coefficient
		outfile << "+1.00000000000e+002 +1.00000000000e+002 +1.00000000000e+002 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PEKS 5\n"; //[MPROP] array containing viscous damping of material for finite elements for each of the properties(for D1PEKS equal 2h E rho, finite element smaller than h is critically damped)
		outfile << "+3.00000000000e+002 +3.00000000000e+002 +3.00000000000e+002 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PELA 5\n";//[MPROP] array containing Lamé elastic constant. (lambda)
		outfile << "+8.65039735565e+009 +8.65039735565e+009 +8.65039735565e+009 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PEMU 5\n";//[MPROP] array containing Lamé elastic constant. (mu)
		outfile << "+1.24481327801e+010 +1.24481327801e+010 +1.24481327801e+010 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PEPE 5\n";//[MPROP] array containing elastic penalty term for contact (usually 2 to 100 times greater than lambda)
		outfile << "+3.00000000000e+011 +3.00000000000e+011 +3.00000000000e+011 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PEPC 5\n";//[MPROP] array containing contact penalty term for contact (usually 2 to 100 times greater than lambda)
		outfile << "+3.00000000000e+011 +3.00000000000e+011 +3.00000000000e+011 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PEPF 5\n"; //[MPROP] array containing static friction coefficient 
		outfile << "+0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PBIF 5\n"; //[MPROP] array containing interface friction coefficient
		outfile << "+6.00000000000e-001 +6.00000000000e-001 +6.00000000000e-001 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PERO 5\n";//[MPROP] array containing density.
		outfile << "+2.70000000000e+003 +2.70000000000e+003 +2.70000000000e+003 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PNAP 5\n";//[MPROP] array containing amplitude of load applied as element surface pressure.
		outfile << "+0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PNAF 5\n";//[MPROP] array containing amplitude factor – all amplitudes are multiplied by this factor.
		outfile << "+0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +1.00000000000e+000\n";
		outfile << "/YD/YDP/D1PNAI 5\n";//[MPROP] array containing the increment of the amplitude factor each time step 
		outfile << "+0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PNAT 5\n";//[MPROP] array containing amplitude of load applied as element surface traction.
		outfile << "+0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PNAX 5\n";//[MPROP] array containing amplitude of force or velocity in local x direction for each node.
		outfile << "+0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PNAY 5\n";//[MPROP] array containing amplitude of force or velocity in local y direction for each node.
		outfile << "+0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PNXX 5\n";//[MPROP] array containing the x component of x axis of local coordinate system for each node.
		outfile << "+1.00000000000e+000 +1.00000000000e+000 +1.00000000000e+000 +1.00000000000e+000 +1.00000000000e+000\n";
		outfile << "/YD/YDP/D1PNXY 5\n";//[MPROP] array containing the y component of x axis of local coordinate system for each node.
		outfile << "+0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PNYX 5\n";//[MPROP] array containing the x component of y axis of local coordinate system for each node.
		outfile << "+0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PNYY 5\n";//[MPROP] array containing the y component of y axis of local coordinate system for each node.
		outfile << "+1.00000000000e+000 +1.00000000000e+000 +1.00000000000e+000 +1.00000000000e+000 +1.00000000000e+000\n";
		outfile << "/YD/YDP/D1PSEM 5\n";//[MPROP] array containing the maximum tensile stretch.
		outfile << "+0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PJRC 5\n";//[MPROP] array containing joint JRC
		outfile << "+1.50000000000e+001 +1.50000000000e+001 +1.50000000000e+001 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PJCS 5\n";//[MPROP] array containing joint JCS
		outfile << "+1.20000000000e+002 +1.20000000000e+002 +1.20000000000e+002 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/D1PJSL 5\n";//[MPROP] array containing joint length
		outfile << "+2.00000000000e-001 +2.00000000000e-001 +2.00000000000e-001 +0.00000000000e+000 +0.00000000000e+000\n";
		outfile << "/YD/YDP/I1PEFR 5\n";//[MPROP] array containing fracture flag. If greater than zero all elements corresponding to this property will be processed as fracturing medium. If it is set to zero for a particular property set, no fracture processing for that property set will take place, i.e. finite elements that are assigned that property will remain non fractured.
		outfile << "0 0 0 0 0\n";
		outfile << "/YD/YDP/I1PEJP 5\n";//[MPROP] array containing joint property(if > 0 joint)
		outfile << "0 2 0 0 0\n";
		outfile << "/YD/YDP/I1PEMB 5\n";//[MPROP] array. If greater than zero boundary nodes are marked, otherwise not.
		outfile << "0 0 0 0 0\n";
		outfile << "/YD/YDP/I1PEMN 5\n";//Numbers of successive mesh refinements.
		outfile << "0 3 0 0 0\n";
		outfile << "/YD/YDP/I1PNFX 5\n";//[MPROP] array. 1 indicates that a force is applied to a node in local x direction, 2 indicates that acceleration is applied to a node and 3 indicates that a velocity is applied to a node. For each node actual force, velocity or acceleration is obtained by multiplying amplitude (such as D1PNAX) with factor (such as D1PNAF).
		outfile << "0 0 0 1 3\n";
		outfile << "/YD/YDP/I1PNFY 5\n";//[MPROP] array. 1 indicates that a force is applied to a node in local y direction, 2 indicates that acceleration is applied to a node and 3 indicates that a velocity is applied to a node. For each node actual force, velocity or acceleration is obtained by multiplying amplitude (such as D1PNAY) with factor (such as D1PNAF).
		outfile << "0 0 0 1 3\n";
		outfile << "/YD/YDP/I1PSDE 5\n";//[MPROP] array containing ID of elastic damage state variable, say zero will indicate that first state variable is elastic damaged, one will indicate that elastic damaged is second state variable, i.e. stored in the second column of the array D2PELST.
		outfile << "0 0 0 0 0\n";
		outfile << "/YD/YDP/I1PTYP 5\n";//[MPROP] array containing the type of each property.Indicates type of object to which this property is associated.For instance : (-1) 2D mechanical x, y d.o.f.node; (1) plane stress elastic triangle; (2) rigid triangle; (3) joint; (4) plain stress softening triangle
		outfile << "1 1 3 -1 -1\n"; 
		outfile.close();
	}
	else
	{
		std::cout << "Unable to open " << User.getfile() << std::endl;
	}
}


void Outfile::writeoutput(UserVariables User, ProblemVariables Problem, BonusVariables Bonus, Material Steel, Material Pvb, Material Glas, Material Rock)
{
	std::ofstream outfile;
	outfile.open(User.getfile(), std::ios::out | std::ios::app | std::ios::binary); ;
	if (outfile.is_open())
	{
		outfile << "/*   Output     */\n";	
		outfile << "/YD/YDO/MOHYS     1\n";//Maximum number of history variables.
		outfile << "/YD/YDO/NOHYS     1\n";//Actual number of history variables.
		outfile << "/YD/YDO/DOHYP +5.000000000000e-004\n";//Output history accuracy.
		outfile << "/YD/YDO/D1OHYC     1\n";//[MOHYS] output history factor to scale time.
		outfile << "+1.000000000000e+000\n";
		outfile << "/YD/YDO/D1OHYF     1\n";//[MOHYS] output history factor to scale state.
		outfile << "+1.000000000000e+000\n";
		outfile << "/YD/YDO/D1OHYS     1\n";//[MOHYS] output history state.
		outfile << "+0.000000000000e+000\n";
		outfile << "/YD/YDO/D1OHYT     1\n";//[MOHYS] output history time.
		outfile << "+0.000000000000e+000\n";
		outfile << "/YD/YDO/D1OHYX     1\n";//[MOHYS] output history x coordinate of the point.
		outfile << "+0.000000000000e+000\n";
		outfile << "/YD/YDO/D1OHYY     1\n";//[MOHYS] output history y coordinate of the point.
		outfile << "+0.000000000000e+000\n";
		outfile << "/YD/YDO/I1OHYT     1\n";//[MOHYS] output history type, i.e. which variable:0 – default (nothing)
			//1 – stress sigmaxx		11 – isotropic elastic damage
			//2 – stress sigmaxy		12 – borehole mass
			//3 – stress sigmayy		13 – borehole pressure
			//7 – velocity				14 – borehole spec.volume
			//8 – velocity x			15 – total kinetic energy
			//9 – velocity y			16 – G2 pressure
			//10 – velocity z			17 – isotropic elastic damage
		outfile << "15\n";
		outfile << "$YDOIT\n"; //It is a command to execute the program.
		outfile << "$YSTOP\n"; //It is a command to stop the program.
		outfile.close();
	}
	else
	{
		std::cout << "Unable to open " << User.getfile() << std::endl;
	}
}

int main(int argc, char **argv)
{
	//arugments from the cmd line: argv[0]=./Ygen, argv[1]=...

	//ready all input files
	std::vector <std::string> geometrydata, userdata, problemdata, bonusdata, steeldata, pvbdata, glasdata, rockdata;
	std::string steel, pvb, glas, rock;
	Inputfile infile;

	geometrydata = infile.readgeometry();
	userdata = infile.readuser();
	problemdata = infile.readproblem();
	bonusdata = infile.readbonus();

	steeldata = infile.readmaterial(steel);
	pvbdata = infile.readmaterial(pvb);
	glasdata = infile.readmaterial(glas);
	rockdata = infile.readmaterial(rock);

	Material Rock(rockdata), Steel(steeldata), Pvb(pvbdata), Glas(glasdata);
	UserVariables User(userdata);
	BonusVariables Bonus(bonusdata, Rock, Steel, Pvb, Glas);
	ProblemVariables Problem(problemdata, Bonus);
	GeometryVariables Geometry(geometrydata);

	Outfile testfile;

	testfile.writecontrol(User, Problem, Bonus, Steel, Pvb, Glas, Rock);
	testfile.writeelements(Geometry, User, Problem, Bonus, Steel, Pvb, Glas, Rock);
	testfile.writeinteractions(User, Problem, Bonus, Steel, Pvb, Glas, Rock);
	testfile.writenodes(User, Problem, Bonus, Steel, Pvb, Glas, Rock);
	testfile.writeproperties(User, Problem, Bonus, Steel, Pvb, Glas, Rock);
	testfile.writeoutput(User, Problem, Bonus, Steel, Pvb, Glas, Rock);
	return 0;
}
