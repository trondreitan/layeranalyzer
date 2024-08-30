#define DETAILED_TIMERS

//#include <lapacke.h>
#include <sys/types.h>

#ifdef MAIN
#include <unistd.h>
#endif // MAIN

#ifndef MAIN // Only for when compiling as R package

#include <Rcpp.h>
#include <RcppCommon.h>

#include <stdlib.h>
#include <stdarg.h>
#include <inttypes.h>

using namespace Rcpp;

//#include <R_ext/Lapack.h> 

#else // MAIN

//#include <R_ext/Lapack.h> 

#endif //MAIN


#include <cmath>

// 
// Made by Trond Reitan, University of Oslo, 15/1-2009.
// LGPL licenced.
//

#include <lapack.h>

// *******************
// *******************
// Definitions:
// *******************
// *******************

#ifndef doubledelete
#define doubledelete(array,len) {if(array) { for(unsigned int count=0;count<(unsigned int)len;count++) if(array[count]) delete [] array[count]; delete [] array; array=NULL;}}
#endif // doubledelete

#ifndef tripledelete
#define tripledelete(array,len,len2) {if(array) { for(unsigned int count=0;count<(unsigned int)len;count++){ if(array[count]){ for(unsigned int count2=0;count2<(unsigned int)len2;count2++) delete [] array[count][count2]; delete [] array[count];}} delete [] array; array=NULL;}}
#endif // tripledelete

#define MAXIM(a, b) ((a) > (b) ? (a) : (b))
#define MINIM(a, b) ((a) < (b) ? (a) : (b))
#define ABSVAL(x) (x>0 ? x : -x)
#define SIGN(x) (x>0 ? 1 : (x==0 ? 0 : -1))

#define MISSING_VALUE -10000000.0  

// two dimenstional structure, used in regression analysis;
struct double_2d
{
  double x,y;
};

// Parameter types: expectancy, characteristic times, diffiusions,
// correlations, external timeseries regression coefficients and
// indicator functions (for grouping sites together, not used in
// the paper). 
enum param_type {MU, LIN_T, DT, SIGMA, CORR, PAIR_CORR, SERIES_CORR, 
		 BETA, INIT, INDICATOR, OBS_SD, TRIG};

// Transformation types for parameters: linear (no transformation),
// logarithmic, logarithm of the inverse, logistic and binary:
enum transform_type {T_LIN, T_LOG, T_LOG_INV, T_INV, T_LOGIST_GLOBAL, 
		     T_LOGIST_PAIR, T_BINARY};


enum INDENTIFICATION_PRIOR_HANDLING {ID_NONE=0,ID_ADD_UPPER=1,ID_SUB_LOWER=2,
				     ID_CUT_LOWER=3,ID_CUT_UPPER=4};
enum REALIZATION_STRATEGY {NO_CENSORING, ASCENDING_CENSORING,
			   DESCENDING_CENSORING};

#define LARGE_ENOUGH_ARRAY 10000

// *******************
// *******************
// Global parameters:
// *******************
// *******************

class series;
class measurement_cluster;
class HydDateTime;

// Parameter info:
char **par_name=NULL; // parameter name
param_type *par_type=NULL; // parameter type
int *par_layer=NULL;  // parameter belonging to which layer
int *par_region=NULL; // parameter belonging to which site
int *par_series=NULL; // parameter belonging to which series

#ifdef DETAILED_TIMERS
long int timers[100][3];
// Timer 1: loglikelihood
// Timer 2: Initialization of loglikelihood
// Timer 3: main loop in loglikelihood
// Timer 4: Kalman smoothing in loglikelihood
// Timer 5: Kalman sampling in loglikelihood
// Timer 9: log_prior
// Timer 10: matrix_inverse. Usses dgesv.
// Timer 11: complex_matrix_inverse
// Timer 12: matrix_eigenvectors. Uses dgeev. Used in Complex_eigenvalues,
//           matrix_determinant, log_matrix_determinant.
// Timer 13: multinormal_sample. Uses dpotrf in get_cholesky.
// Timer 14: pdf_multinormal, pdf_multi_t (all)
// Timer 15: mean_of_vectors/estimated_variance
// Timer 20: Initialization of MCMC:
// Timer 21: Burnin
// Timer 22: MCMC sampling
// Timer 23: Bayesian model log-likelihood or ML optimization
// Timer 24: test_OU
// Timer 25: runs_test
// Timer 26: analyze_residuals
#endif // DETAILED_TIMERS



// Transformation of parameters:

// transformation type of each parameter contained here:
transform_type *par_trans_type=NULL;

REALIZATION_STRATEGY real_strat=NO_CENSORING;

series *ser=NULL;
unsigned int num_series=0;

unsigned int num_smooth=10; // Smoothings per iteration

INDENTIFICATION_PRIOR_HANDLING id_strategy=ID_NONE;

unsigned int numsites; // Number of sites
unsigned int num_tot_layers; // total number of layers over all series
bool some_pairwise_correlations=false;

double p_pos_site_sigma2=0.0;
double p_pos_series_sigma2=0.0;

// Extra-layer instantaneous correlations specifications:
unsigned int num_series_corr=0;
int *corr_from_series=NULL, *corr_to_series=NULL, 
  *corr_from_layer=NULL, *corr_to_layer=NULL;
int *corr_from_index=NULL, *corr_to_index=NULL; // series+layer indexes
double *series_corr=NULL;


int use_indicator=0;
int *indicator_array=NULL;
int useext=0; // use external times series in the main series?

measurement_cluster *meas_tot=NULL, *meas_smooth=NULL;
unsigned int meas_tot_len=0, meas_smooth_len=0;

double *obs_corr_t=NULL, **obs_corr=NULL;
int *obs_corr_site=NULL, obs_corr_len=0;

unsigned int num_states=1;
unsigned int state_series[LARGE_ENOUGH_ARRAY]; // which serie does a state variable belong to?
unsigned int state_layer[LARGE_ENOUGH_ARRAY]; // which layer does a state variable belong to?
unsigned int series_state_start[100]; // which state index does this series start with?

unsigned int numpar=0; // number of parameters


unsigned int len_x_k_s_kept=0, size_x_k_s_kept=0;
double **x_k_s_kept=NULL, ***P_k_s_kept=NULL;
double *t_k_smooth=NULL;
HydDateTime *dt_k_smooth=NULL;

double **x_k_realized=NULL;
double ***x_k_realized_all=NULL;
unsigned int numit_realization=100, numit_realizations_made=0;

bool nodata=false;

// Extra-layer causal specifications:
unsigned int num_series_feed=0;
int *feed_from_series=NULL, *feed_to_series=NULL, 
  *feed_from_layer=NULL, *feed_to_layer=NULL;
int *feed_symmetric=NULL;
double *beta_feed=NULL;
bool is_complex=false;
double cycle_time[10];

// Measurement data and external data series:
unsigned int ext_len=0;
double_2d *extdata=NULL;


// Run options:
int silent=0; // indicates silent modues. Default switched off.
int talkative_likelihood=0;
int talkative_burnin=0;


// *************************************************
// Purging mechanism for global variables, so that
// they can be reused.
// *************************************************

void reset_global_variables(void);


// ****************************************************
// General purpose date-time format
// Taken from LGPL-licensed hydrasub-library, files 
// ctime, date and date_time. Developed
// at NVE 1996-2016 (http://www.nve.no).
// ****************************************************

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cstdio>


//#ifdef MAIN
using namespace std;
//#endif // MAIN

typedef unsigned short hourType;
typedef unsigned short minuteType;
typedef unsigned short secondType;
typedef unsigned int totalSecondsType;

int legalTime(const short h, const short m, const short s);

void skipDelim(istream &is) 
{
  // Help function to read and discard the delimiting character
  char c;
  if (!is.good()) return;
  is >> c;
  while (is.good() && !isalnum(c)) is >> c;
  if (is.good()) is.putback(c);
}

class CTime 
{    
protected:
    totalSecondsType tsec;
    CTime(totalSecondsType t) : tsec(t) {};
    
public:
    void totSec2time(hourType&, minuteType&, secondType&) const;
    
    /* Convert seconds to hhmmss CTime (hour:minute:sec). This function is 
       not valid for CTime whose total seconds is greter 
       than 2**sizeof(int)-1. */
    
    totalSecondsType time2totSec(hourType, minuteType, secondType) const;
    
    /* Convert an hhmmss CTime (hour:minute:sec) to seconds. 
       This function is not valid for CTime whose total seconds is greter 
       than 2**sizeof(int)-1. */
    
public:
    CTime () {tsec=0;};
    CTime (istream&);
    CTime (const CTime &);
    CTime (const char *);
    CTime (hourType, minuteType, secondType=0);

    void now();
    int legal() const;

    // assignment operator
    CTime& operator= (const CTime &t2) { 
	tsec = t2.tsec; 
	return *this;
    };
    
    // assignment operators
    friend int operator== (const CTime &t1, const CTime &t2) { 
	return int(t1.tsec == t2.tsec); };           
    friend int operator!= (const CTime &t1, const CTime &t2) { 
	return int(t1.tsec != t2.tsec); };
    friend int operator<(const CTime &t1, const CTime &t2)  { 
	return int(t1.tsec < t2.tsec); };
    friend int operator<=(const CTime &t1, const CTime &t2) { 
	return int(t1.tsec <= t2.tsec); };
    friend int operator> (const CTime &t1, const CTime &t2)  { 
	return t2 < t1; };
    friend int operator>= (const CTime &t1, const CTime &t2) { 
	return t2 <= t1; };
    
    friend CTime operator+ (const CTime &t, int sec) 
	{ return CTime(t.tsec + sec); };
    friend CTime operator+ (const CTime &t1, const CTime &t2)
	{ return CTime(t1.tsec + t2.tsec); };
    friend void operator+= (CTime &t, int sec)  
	{ t.tsec += sec; };
    
    friend int operator- (const CTime &t1, const CTime &t2)
	{ return (int) (t1.tsec - t2.tsec); };
    friend CTime operator- (CTime &t, int sec) 
	{ return CTime(t.tsec - sec); };
    friend void operator-= (CTime &t, int sec)
	{ t.tsec = t.tsec -sec; };
    // friend void operator-= (CTime &t, const CTime &t2)
    //    { t.tsec = t.tsec -t2.tsec; };
    
    friend int between(const CTime &t3, const CTime &t1, const CTime &t2)
	{ return int(t3.tsec >= t1.tsec && t3.tsec <= t2.tsec); };
    
    friend ostream& operator<<(ostream&, const CTime &);
    //void print();
    
    hourType getHour() const;
    minuteType getMinute() const;
    secondType getSecond() const;
    char *getChTime() const;

    int minofday() const;		// [0 - 1399]

    int secofday() const;		// [0 - 86399]

      // help function for CTime(istream&)
    totalSecondsType parseTime(istream&)const;
};




class HydDateTime;

typedef unsigned short dayType;
typedef unsigned short monthType;
typedef unsigned short yearType;
typedef unsigned int   julType;
enum WEEKDAY {MONDAY=0, TUESDAY=1, WEDNESDAY=2, THURSDAY=3, 
	      FRIDAY=4, SATURDAY=5, SUNDAY=6};

int leap(yearType y);
int legalHydDate(dayType dd, monthType mm, yearType yyyy);


class HydDate 
{  
  friend class HydDateTime;
protected:
  julType julnum;
  HydDate(julType j) { julnum = j; };
    
public:
  julType date2jday(void);
  julType date2jday(dayType, monthType, yearType) const;
  /* Convert a Gregorian Calender date to its corresponding 
     Julian day number. Gregorian Calender started on Sep. 14, 1752 
     (we use this to denote an illegal date).
     This function is not valid before that date. */
  
  void jday2date(dayType&, monthType&, yearType&) const;
  /* Convert a Julian day number to its corresponding Gregorian
     Calender date. Gregorian Calender started on Sep. 14, 1752.
     This function is not valid before that date. */
  
public:
  HydDate();
  HydDate(istream& is);
  HydDate(const char *);
  HydDate(const HydDate &);
  HydDate(yearType y, monthType m, dayType d = 1); 
  
  void now(); // sets the date to today
  int legal() const; /* returns true if a legal date.
			OBS: illegal date is stored as 1752.sep.14*/
  
  // assignment operator
  HydDate& operator=(const HydDate &d1) {
    julnum = d1.julnum; return *this;};
  HydDate& operator=(const HydDateTime &d1);
  
  // logical operators
  friend int operator== (const HydDate &d1, const HydDate &d2) { 
    return int(d1.julnum == d2.julnum); };           
  friend int operator!= (const HydDate &d1, const HydDate &d2) { 
    return int(d1.julnum != d2.julnum); };
  friend int operator<(const HydDate &d1, const HydDate &d2)  { 
    return int(d1.julnum < d2.julnum); };
  friend int operator<=(const HydDate &d1, const HydDate &d2) { 
    return int(d1.julnum <= d2.julnum); };
  friend int operator> (const HydDate &d1, const HydDate &d2)  { 
    return int(d2 < d1); };
  friend int operator>= (const HydDate &d1, const HydDate &d2) { 
    return int(d2 <= d1);  };
    
  // arithmetical operators
    
  // adds some days 'dd' to a date and returns the resulting date
  friend HydDate operator+ (const HydDate &d, int dd) 
    { return HydDate(d.julnum + dd); };
  friend HydDate operator+ (int dd, const HydDate &d)
    { return HydDate(dd + d.julnum); };
  
  //adds days 'dd' to the given date 'd';
  //note: if you don't want the original date changed, use + mentioned above
  friend void operator+= (HydDate &d, int dd) { d.julnum += dd; };
    
  //returns the difference between 2 dates in number of days
  friend int operator- (const HydDate &d1, const HydDate &d2) 
    { return (d1.julnum - d2.julnum); };
    
  // substract 'dd' days from a date 'd' and returns the resulting date
  friend HydDate operator- (const HydDate &d, int dd) 
    { return HydDate(d.julnum - dd); };
    
  friend void operator-= (HydDate &d, int dd) { d.julnum -= dd; };
    
  // others
  friend int between(const HydDate &d3, const HydDate &d1,  const HydDate &d2)
    {return int(d3.julnum >= d1.julnum && d3.julnum <= d2.julnum); };
    
  friend int overlap_date(const HydDate &start1, const HydDate &end1, 
                          const HydDate &start2, const HydDate &end2);

  friend ostream& operator<<(ostream&, HydDate &);
    //  void print();
    
  char *charHydDate() const;  // yyyy.mm.dd
  yearType getYear() const;
  monthType getMonth() const;
  dayType getDay() const;
  int getDayOfWeek() const; // sunday = 0, saturday = 6
    
  int leap() const; // whether the current date-object's year is leap
    
  julType parseHydDate(istream&) const; // help function for HydDate(istream&)
    
  int dayOfYear() const;    // returns the corresponding day of the year
			    // starting at 1

					     
  int dayIndex() const;	  // returns the corresponding day of the year
			  // starting at 1, ingoring leap days


  int daysInYear() const; // returns no. of days in the year: 365 or 366
  int firstDayOfMonth() const;	// returns the first day of the month
  int noDaysInMonth() const;	// returns no. of days in the month
    
  int idxOfYear() const;		// return pos in year starting at 

  int Get_weekday_index(void);
  WEEKDAY Get_weekday(void);
};

HydDate NoHydDate;


int legalHydDateTime(yearType y, monthType m, dayType d, short h, short min, 
		  short s);

class HydDateTime : public  HydDate, public CTime {


    HydDateTime (julType j, totalSecondsType ts) : HydDate(j), CTime(ts) { };    
    
public:
    HydDateTime () ;
    HydDateTime (const char *date, const char* time);
    HydDateTime (yearType yy, monthType mm = 1, dayType dd = 1, hourType hh = 0, 
	      minuteType mn = 0, secondType second = 0);

    HydDateTime (const HydDate &d);
    HydDateTime(const HydDateTime& d);
    HydDateTime(const char *datetime);

    HydDateTime(char *datetime,int form); // se syCh(int form)

    void now();				// sets the datetime to now.
    int legal() const;		// checks for legal datetime.

    HydDateTime Floor(int precision_min); // round down
    // to the given precision (in minutes)
    HydDateTime Ceiling(int precision_min); // round up
    // to the given precision (in minutes)

    // difference in minutes between this datetime and the second, 
    // given datetime
    int difference_minutes(HydDateTime &dt2);	 

    HydDateTime &add_minutes(int num_minutes);

    // assignment operator
    HydDateTime& operator=(const HydDateTime &dt2);
    HydDateTime& operator=(const HydDate &dt2);
    
    // logical operators
    friend int operator== (const HydDateTime &dt1, const HydDateTime &dt2);
    friend int operator!= (const HydDateTime &dt1, const HydDateTime &dt2);
    friend int operator<  (const HydDateTime &dt1, const HydDateTime &dt2);
    friend int operator<= (const HydDateTime &dt1, const HydDateTime &dt2);
    friend int operator>  (const HydDateTime &dt1, const HydDateTime &dt2);
    friend int operator>= (const HydDateTime &dt1, const HydDateTime &dt2);
    
    // arithmetical operators
    
    friend HydDateTime operator+ (const HydDateTime &dt, int sec);
    friend HydDateTime operator+ (int sec, const HydDateTime &dt);
    
    friend unsigned int operator- (const HydDateTime &dt1, const HydDateTime &dt2);
    // this func. returns the absolute difference between the two times
    // in seconds
    
    friend HydDateTime operator- (const HydDateTime &d, int sec);
    friend HydDateTime operator- (int sec, const HydDateTime &d);
    // this functions substract 'sec' seconds from the datetime 'd'
    // and returns the substracted datetime
    
    friend void operator-= (HydDateTime &dt, int sec);
    friend void operator+= (HydDateTime &dt, int sec);
    
    // others
    
    // d1 <= dt2 && dt2 <= dt3
    friend int between(const HydDateTime &d1, const HydDateTime &dt2,  
		       const HydDateTime &dt3);
    
    // Does two periods overlap?
    friend int overlap(const HydDateTime &start1, const HydDateTime &end1,  
		       const HydDateTime &start2, const HydDateTime &end2);

    friend ostream& operator<<(ostream&, const HydDateTime &);
    void Print();
    void getHydDateTime(yearType &yyyy , monthType &mm, dayType &dd, 
		     hourType &hh, minuteType &mn, secondType &ss);

    char *syCh(int form = 0) const;		// Sybase format
    // 0 - MM/dd/yyyy hh:mm:ss
    // 1 - dd.MM.yyyy hh:mm
    // 2 - yyyyMMddhhmm
    // 3 - yyyy.MM.dd hh:mm
    // 4 - dd.MM.yyyy
    // 5 - yyyyMMdd/hhmm
    // 6 - dd.MM hh:mm
    // 7 - dd.mm.yyyy hh:mm:ss
    // Changed '/' and '-' to '.' 23/1-2013 by Trond Reitan

    
    

    int  isNull() const;


    HydDateTime StartOfDay() const;		// yyyymmdd 00:00:00
    HydDateTime EndOfDay() const;		        // yyyymmdd 23:59:00
    
    HydDateTime StartOfMonth() const;		// yyyymmdd 00:00:00  	   
    HydDateTime EndOfMonth() const;	        // yyyymmdd 23:59:0 	   

    HydDateTime StartOfYear() const;		// yyyymmdd 00:00:00
    HydDateTime EndOfYear() const;		        // yyyymmdd 23:59:00			     
    double as_floating_point_year(void);
						     
};


static HydDateTime NoHydDateTime(1753, 1, 1, 0, 0, 0); //erikt


int operator== (const HydDateTime &dt1, const HydDateTime &dt2) 
{ 
  return int(dt1.julnum == dt2.julnum && dt1.tsec == dt2.tsec); 
}

int operator!= (const HydDateTime &dt1, const HydDateTime &dt2) 
{ 
  return int(dt1.julnum != dt2.julnum || dt1.tsec != dt2.tsec); 
}


int operator> (const HydDateTime &dt1, const HydDateTime &dt2)  
{ 
  return dt2 < dt1; 
}

int operator>= (const HydDateTime &dt1, const HydDateTime &dt2) 
{ 
  return dt2 <= dt1;  
}

HydDateTime maxdt(HydDateTime dt1,HydDateTime dt2)
{
  if(dt1!=NoHydDateTime && (dt1>dt2 || dt2==NoHydDateTime))
    return dt1;
   else
     return dt2;
}

HydDateTime mindt(HydDateTime dt1,HydDateTime dt2)
{
   if(dt1!=NoHydDateTime && (dt1<dt2 || dt2==NoHydDateTime))
     return dt1;
   else
     return dt2;
}

int between(const HydDateTime &dt1, const HydDateTime &dt2,  
		       const HydDateTime &dt3) 
{
  return int( (dt2!=NoHydDateTime) && ( (HydDate) dt2!=NoHydDate) && 
	      (dt1==NoHydDateTime || dt1 <= dt2) && 
	      (dt3==NoHydDateTime || dt2 <= dt3) ); 
}


int legalTime(const short h, const short m, const short s) {
    return (h<0 || h>23 || m < 0 || m > 59 || s < 0 || s > 59)?
      0 : 1;
}

totalSecondsType CTime::time2totSec(hourType h, minuteType m, secondType s) const {
    /* Convert an hhmmss time (hour:minute:sec) to seconds. 
       This function is not valid for time whose total seconds is greater 
       than 2**sizeof(unsigned int). */
    if (!legalTime(h,m,s)) {
	h = 24; m = 60; s = 60;
    };
    return (s + m*60 + h*60*60);
    
} /* time2totSec */


void CTime::totSec2time(hourType &h, minuteType &m, secondType &s) const {
    /* Convert seconds to hhmmss time (hour:minute:sec). This function is 
       not valid for time whose total seconds is greter 
       than 2**sizeof(unsigned int). */
    
    h = tsec/3600;
    m = (tsec % 3600) / 60;
    s = (tsec % 3600) % 60;
    
} /* totSec2time */

void CTime::now() {
    time_t clk = time(0);
    tm* now = NULL;
    now = localtime(&clk);
    tsec = time2totSec(now->tm_hour, now->tm_min, now->tm_sec);
}

CTime::CTime(const CTime &t) {
    tsec = t.tsec;
}


CTime::CTime(hourType h, minuteType m, secondType s) {
    tsec = time2totSec(h,m,s);
}

CTime::CTime(const char *time) {
    static char tem[3];
    static hourType hh = 0; static minuteType mm = 0; static secondType ss = 0;
    tem[2] = '\0';
    hh = (hourType)atoi(strncpy(tem, time, 2));
    mm = (minuteType)atoi(strncpy(tem, time + 2, 2));
    ss = 0;
    if (strlen(time) > 4)
      ss = atoi(strncpy(tem, time + 4, 2));
    if(legalTime(hh,mm,ss))
      tsec = time2totSec(hh, mm, ss);
    else
      tsec = time2totSec(24, 60, 60);
}

int CTime::legal() const {
    static hourType hh = 0;  static minuteType mm = 0;  static secondType ss = 0;
    totSec2time(hh, mm, ss);
    return (::legalTime(hh,mm,ss));
}

totalSecondsType CTime::parseTime(istream& is) const {
    /* parse time of the format hh:mm:ss
       Delimiter need not be ':',
       it can be any of the alphanumeric characters */
    
    unsigned h = 0, m = 0, s = 0;
    
    if (is.good()) {
	skipDelim(is);
	is >> h;
	skipDelim(is);
	if (is.eof()) return 0;
	
	is >> m;
	skipDelim(is);
	if (is.eof()) return 0;
	is >> s;
    }
    if (!is.good()) return 0;
    return CTime((hourType)h,m,s).tsec;
}

// end parseTime

CTime::CTime(istream &is) { 
    tsec = parseTime(is); 
}


ostream& operator<<(ostream& os, const CTime &t) 
{    
  static hourType hh = 0;  
  static minuteType mm = 0;  
  static secondType ss = 0;
  static char buffer[100];
  
  t.totSec2time(hh, mm, ss);
  snprintf(buffer, 99, "%02d:%02d", hh, mm);
  os << buffer;
  
  return os;
}


char *CTime::getChTime() const
{
  static char ctime[100];
  static hourType h = 0;  static minuteType mn = 0;  static secondType s = 0;
  
  totSec2time(h, mn, s);
  snprintf(ctime, 99, "%02d:%02d", h, mn);
  return ctime;
 }

hourType CTime::getHour() const {
    static hourType hh = 0;  static minuteType mm = 0;  static secondType ss = 0;
    totSec2time(hh, mm, ss);
    return hh;
}

minuteType CTime::getMinute() const {
    static hourType hh = 0;  static minuteType mm = 0;  static secondType ss = 0;
    totSec2time(hh, mm, ss);
    return mm;
}

secondType CTime::getSecond() const {
    static hourType hh = 0;  static minuteType mm = 0;  static secondType ss = 0;
    totSec2time(hh, mm, ss);
    return ss;
}


int CTime::minofday() const {
    return tsec / 60;
}

int CTime::secofday() const {
    return tsec;
}


// C and Fortran interface:
extern "C" {

int monthsize(char *datestr)
{
  char str[20];
  
  if(strlen(datestr)<6)
    return 0;
  
  strncpy(str, datestr, 6);
  strcpy(str+6,"01");
  
  HydDate dt(str);
  
  if(!dt.legal())
    return 0;
  
  int size=(int) dt.noDaysInMonth();
  
  return size;
}

int monthsize_(char *datestr)
{
  return monthsize(datestr);
}

} // extern "C"

/* day_tab includes the days of the months in year/leap year */
static int day_tab[2][13] = {
    {0,31,28,31,30,31,30,31,31,30,31,30,31},
    {0,31,29,31,30,31,30,31,31,30,31,30,31}
};

static int firstdayinmonthtab[2][13] = {
    {0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334},
    {0, 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335}
};


int leap(yearType y) 
{
  return  int((y%4 == 0 && y%100 != 0) || (y%400 == 0 && y%4000!=0));
}


int legalHydDate(dayType dd, monthType mm, yearType yyyy) 
{
  // illegal month or date < 1
  if (mm < 1 || mm > 12 || dd < 1 || dd > 31 || yyyy < 1753) 
    return 0; 
  
  if (dd > day_tab[leap(yyyy)][mm])
    return 0; // illegal date
  return 1;
}


julType HydDate::date2jday(void)
{
  return julnum;
}


// The algorithm used in the following 2 functions was described in 
// Communications of the ACM, Vol. 6, No. 8 (Aug.63), p.444.

julType HydDate::date2jday(dayType d, monthType m, yearType y) const 
{
  /* Convert a Gregorian Calender date to its corresponding 
     Julian day number. Gregorian Calender started on  14,Sep. 1752.
     This function is not valid before that date. 
     NB: WE USE 1.JAN.1753 TO DENOTE ILLEGAL DATE.
     FOR EXAMPLE, IF AN ILLEGAL DATE, 35 JAN 1994
     IS GIVEN TO CREATE A DATE VARIABLE, THEN THE
     VALUE STORED IN THERE IS 1.JAN.1753. */
  
  if (!legalHydDate(d,m,y)) 
    {
      d = 1; m = 1; y = 1753;
    }

  static unsigned int c = 0, ya = 0;
  if (m>2) 
    m -= 3;
  else 
    {
      m += 9;
      y--;
    }
  c = y / 100;
  ya = y - 100*c;
  return ((146097*c)>>2) + ((1461*ya)>>2) + (153*m +2)/5 + d + 1721119;
} /* date2jday */


void HydDate::jday2date(dayType &d, monthType &m, yearType &y) const 
{
  /* Convert a Julian day number to its corresponding Gregorian
     Calender date. Gregorian Calender started on Sep. 14, 1752.
     This function is not valid before that date. */
  
  julType j = julnum - 1721119;
  y = (yearType) (((j<<2) - 1) / 146097);
  j = (j<<2) - 1 - 146097*y;
  d = (dayType) (j>>2);
  j = ((d<<2) + 3) / 1461;
  d = (dayType) ((d<<2) + 3 - 1461*j);
  d = (d+4)>>2;
  m = (5*d -3)/153;
  d = 5*d -3 - 153*m;
  d = (d+5)/5;
  y = (yearType)(100*y + j);
  
  if (m < 10)  
    m += 3;
  else 
    {
      m -= 9;
      y++;
    }
} /* jday2date */


HydDate::HydDate() 
{
  julnum = 2361331;		// Set to 1 JAN 14, which is NoHydDate
}


HydDate::HydDate(const HydDate &d) 
{
  julnum = d.julnum;
}


HydDate::HydDate(yearType y, monthType m, dayType d) 
{
  julnum = date2jday(d,m,y);
}


HydDate::HydDate(const char *date) 
{	
  //yyyymmdd
  if (strcmp(date, "null") == 0) 
    julnum = date2jday(1, 1, 1753);
  else 
    {
      char temp[5]; char temp1[3];
      static dayType dd = 0;  static monthType mm = 0;  
      static yearType yyyy = 0;
      
      temp[4] = temp1[2] = '\0';
      yyyy =  (yearType) atoi(strncpy(temp,date,4));
      mm = (monthType) atoi(strncpy(temp1,date+4,2));
      dd = (dayType) atoi(strncpy(temp1,date+6,2));
      
      julnum = date2jday(dd, mm, yyyy);
    }
}



julType HydDate::parseHydDate(istream& is) const 
{
  /* parse dates of the format yyyy.mm.dd
     or dd.mm.yyyy. Delimiter need not be '.',
     it can be any of the alphanumeric character */
  static dayType d = 0;  static monthType m = 0;  static yearType y = 0;  
  static unsigned temp = 0;
    
  if (is.good()) 
    {
      skipDelim(is);
      is >> temp;
      skipDelim(is);
      if (is.eof()) return 0;
      if (temp / 100 > 0) 
	{   
	  // date formate is: yyyy.mm.dd
	  y = temp;
	  is >> m;
	  skipDelim(is);
	  if (is.eof()) return 0;
	  is >> d;
	}
      else 
	{
	  // date formate is: dd.mm.yyyy
	  d = temp;
	  is >> m;
	  skipDelim(is);
	  if (is.eof()) return 0;
	  is >> y;
	}
    }
  if (!is.good()) 
    return 0;
  return HydDate((yearType)y,m,d).julnum;
} // end parseHydDate


HydDate::HydDate(istream &is) 
{   
  julnum = parseHydDate(is); 
}


ostream& operator<<(ostream& os,  HydDate &d) 
{
  if(d == NoHydDate) 
    os << "NULL";
  else 
    {
      static dayType dd = 0; 
      static monthType mm = 0; 
      static yearType yyyy = 0;
      static char buffer[100];

      d.jday2date(dd, mm, yyyy);
      snprintf(buffer, 99, "%02d.%02d.%d ", dd, mm, yyyy);
      os << buffer;
    }
  
  return os;
}


int HydDate::Get_weekday_index(void)
{
  static HydDate monday(1999,1,25);

  int weekday, diff = (int) abs((int)(julnum - monday.julnum));
  if(julnum < monday.julnum)
    weekday = (7 - diff % 7) % 7;
  else
    weekday = diff % 7;
  return weekday;
}


WEEKDAY HydDate::Get_weekday(void)
{
  return (WEEKDAY) Get_weekday_index();
} 



void HydDate::now() 
{ 
  time_t clk = time(0);
  tm* now;
  now = localtime(&clk);
  if(now->tm_year<50)
    julnum = date2jday(now->tm_mday, now->tm_mon+1, now->tm_year + 2000);
  else
    julnum = date2jday(now->tm_mday, now->tm_mon+1, now->tm_year + 1900);
}


int HydDate::legal() const 
{
  return (*this != NoHydDate) ? 1 : 0;
}


HydDate& HydDate::operator=(const HydDateTime &d1) 
{
  julnum = d1.julnum; 
  return *this;
}


yearType HydDate::getYear() const 
{
  static dayType dd = 0; static monthType mm = 0; static yearType yyyy = 0;
  jday2date(dd, mm, yyyy);
  return yyyy;
}


monthType HydDate::getMonth() const 
{
  static dayType dd = 0;  static monthType mm = 0;  static yearType yyyy = 0;
  jday2date(dd, mm, yyyy);
  return mm;
}


dayType HydDate::getDay() const 
{
  static dayType dd = 0;  static monthType mm = 0;  static yearType yyyy = 0;
  jday2date(dd, mm, yyyy);
  return dd;
}

int HydDate::leap() const 
{
  return (::leap(getYear())) ? 1 : 0;
}






int HydDate::dayOfYear() const 
{
  static dayType dd = 0; static monthType mm = 0; static yearType yy = 0;
  jday2date(dd, mm, yy);
  
  return firstdayinmonthtab[::leap(yy)][mm] + dd;
}


int HydDate::dayIndex() const
{
  static dayType dd = 0; static monthType mm = 0; static yearType yy = 0;
  jday2date(dd, mm, yy);
  
  if(mm==2 && dd==29)
    return -1;
  
  return firstdayinmonthtab[0][mm] + dd;
}

int HydDate::daysInYear() const 
{
  return  (leap()) ? 366 : 365;
}


int HydDate::firstDayOfMonth() const 
{
  static dayType dd = 0; static monthType mm = 0; static yearType yyyy = 0; 
  int d = 0, i;
  jday2date(dd, mm, yyyy);
  for (i = 1; i < mm; i++) 
    d += day_tab[leap()][i];
  d++;
  return d;
}


int HydDate::noDaysInMonth() const 
{
  static dayType dd = 0; static monthType mm = 0; static yearType yyyy = 0;
  jday2date(dd, mm, yyyy);
  return day_tab[leap()][mm];
}



char *HydDate::charHydDate() const 
{
  char *pt = new char[100];
  if(*this == NoHydDate)
    snprintf(pt, 99, "null");
  else
    snprintf(pt, 99, "%d.%d.%d", getYear(), getDay(), getMonth());
  return (pt);
}


int HydDate::idxOfYear() const 
{
  static dayType dd = 0; static monthType mm = 0; static yearType yy = 0;
  jday2date(dd, mm, yy);
  
  return firstdayinmonthtab[::leap(yy)][mm] + dd - 1; 
}


int HydDate::getDayOfWeek() const 
{
  // 15. sept. 1752 was a Monday (day 1), so:
  return (julnum - date2jday(15,9,1752) +1) %7;
}
	

int overlap_date(const HydDate &start1, const HydDate &end1, 
		 const HydDate &start2, const HydDate &end2)
{
  if(start1!=NoHydDate && end2!=NoHydDate && start1>end2)
    return 0;
  
  if(end1!=NoHydDate && start2!=NoHydDate && end1<start2)
    return 0;

  return 1;
}



const char* Month[12] = 
{
  "Januar", "Februar", "Mars", "April", "Mai", "Juni", "Juli", 
  "August","September", "Oktober", "November", "Desember"
};

const char* Month_eng[12] = 
{
  "January", "February", "March", "April", "May", "June", "July", 
  "August","September", "October", "November", "December"
};

const char* Month_short[12] = 
{
  "Jan", "Feb", "Mar", "Apr", "Mai", "Jun", "Jul", 
  "Aug","Sep", "Okt", "Nov", "Des"
};

const char* Month_eng_short[12] = 
{
  "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", 
  "Aug","Sep", "Oct", "Nov", "Dec"
};

int legalHydDateTime(yearType y, monthType m, dayType d,
		  short h, short min, short s) {
  return (legalHydDate(y,m,d) && legalTime(h,min,s)) ? 
    1 : 0;
}

HydDateTime::HydDateTime () : HydDate(), CTime() 
{
  // NOP
} 

HydDateTime::HydDateTime (const char *date, const char* time) : HydDate(date), 
  CTime(time)
{
  // NOP
}

HydDateTime::HydDateTime (yearType yy, monthType mm, dayType dd, 
		    hourType hh, minuteType mn, secondType ss) : 
    HydDate(yy, mm, dd), CTime(hh, mn, ss) 
{
  // NOP
}

HydDateTime::HydDateTime (const HydDate &d) : HydDate(d), CTime((hourType)0, 0, 0)
{
  // NOP
}

HydDateTime::HydDateTime(const HydDateTime& d) :  HydDate(d.julnum), CTime(d.tsec)
{
  // NOP
}

// 199401011200
// 19940101120000
HydDateTime::HydDateTime(const char *datetime)  
{    
  char temp[5]; temp[4] = '\0'; 
  char temp1[3]; temp1[2] = '\0';
  
  yearType  y = atoi(strncpy(temp, datetime, 4));
  monthType m = atoi(strncpy(temp1, datetime+4, 2));
  dayType   d = atoi(strncpy(temp1, datetime+6, 2));
  HydDate::julnum = date2jday(d, m, y);
  
  short h, mi, sec = 0;
  h = atoi(strncpy(temp1, datetime+8, 2));
  mi = atoi(strncpy(temp1, datetime+10, 2));
  if (strlen(datetime) > 12)
    sec = atoi(strncpy(temp1, datetime+12, 2));
  CTime::tsec = time2totSec(h, mi, sec);
    
}

char *nextdigits(char *str,int num)
{
  int i=0,onlydigits=1;
  int strmax=strlen(str);

  while(((i<num && onlydigits) || !isdigit(str[i])) && i<strmax)
    {
      i++;
      if(!isdigit(str[i]))
	onlydigits=0;
    }
  
  if(!isdigit(str[i]))
    return NULL;
  return str+i;
}

HydDateTime::HydDateTime(char *datetime,int form)  
{    
  char temp[5]; temp[4] = '\0'; 
  char temp1[3]; temp1[2] = '\0';
  char *ptr1,*ptr2,*ptr3,*ptr4,*ptr5;
  yearType  y=2000;
  monthType m=1;
  dayType d=1;
  short h=0,mi=0,sec=0;

  if(!datetime)
    return;

  switch(form)
    {
    case 0:
      ptr1=strstr(datetime,"/")+1;
      if(!ptr1)
	{
	  y=0;m=2;d=30;
	  return;
	}
      ptr2=strstr(ptr1,"/")+1;
      if(!ptr2)
	{
	  y=0;m=2;d=30;
	  return;
	}
      sscanf(datetime,"%hu",&m);
      sscanf(ptr1,"%hu",&d);
      sscanf(ptr2,"%hu",&y);
      ptr3=strstr(ptr2," ")+1;
      if(ptr3)
	{
	  sscanf(ptr3,"%hd",&h);
	  ptr4=strstr(ptr3,":")+1;
	  if(ptr4)
	    {
	      sscanf(ptr4,"%hd",&mi);
	      ptr5=strstr(ptr4,":")+1;
	      if(ptr5)
		sscanf(ptr5,"%hd",&sec);
	    }
	}
      break;
    case 1:
    case 7:
      d=atoi(strncpy(temp1, datetime, 2));
      ptr1=nextdigits(datetime,2);
      if(!ptr1)
	{
	  y=0;m=2;d=30;
	  return;
	}
      m=atoi(strncpy(temp1, ptr1, 2));
      ptr2=nextdigits(ptr1,2);
      if(!ptr2)
	{
	  y=0;m=2;d=30;
	  return;
	}
      y=atoi(strncpy(temp, ptr2, 4));
      ptr3=nextdigits(ptr2,4);
      if(!ptr3)
	break;
      h=atoi(strncpy(temp1, ptr3, 2));
      ptr4=nextdigits(ptr3,2);
      if(!ptr4)
	break;
      mi=atoi(strncpy(temp1, ptr4, 2));
      ptr5=nextdigits(ptr4,2);
      if(!ptr5)
	break;
      sec=atoi(strncpy(temp1, ptr5, 2));
      /*
      y = atoi(strncpy(temp, datetime+6, 4));
      m = atoi(strncpy(temp1, datetime+3, 2));
      d = atoi(strncpy(temp1, datetime, 2));
      h = atoi(strncpy(temp1, datetime+11, 2));
      mi = atoi(strncpy(temp1, datetime+14, 2));
      sec = 0;
      */
      break;
    case 2:
      y = atoi(strncpy(temp, datetime, 4));
      m = atoi(strncpy(temp1, datetime+4, 2));
      d = atoi(strncpy(temp1, datetime+6, 2));
      h = atoi(strncpy(temp1, datetime+8, 2));
      mi = atoi(strncpy(temp1, datetime+10, 2));
      sec = 0;
      break;
    case 3:
      y = atoi(strncpy(temp, datetime, 4));
      m = atoi(strncpy(temp1, datetime+5, 2));
      d = atoi(strncpy(temp1, datetime+8, 2));
      h = atoi(strncpy(temp1, datetime+11, 2));
      mi = atoi(strncpy(temp1, datetime+14, 2));
      sec = 0;
      break;
    case 4:
      strncpy(temp, datetime+6, 4);
      temp[4]='\0';
      y = atoi(temp);
      strncpy(temp1, datetime+3, 2);
      temp1[2]='\0';
      m = atoi(temp1);
      strncpy(temp1, datetime, 2);
      temp1[2]='\0';
      d = atoi(temp1);
      h = 0;
      mi = 0;
      sec = 0;
      break;
    }

  HydDate::julnum = date2jday(d, m, y);
  CTime::tsec = time2totSec(h, mi, sec);
}



 
void HydDateTime::now() 
{
  HydDate::now();
  CTime::now();
}

int HydDateTime::legal() const 
{
  return (HydDate::legal() && CTime::legal()) ? 1: 0;
}

HydDateTime& HydDateTime::operator=(const HydDateTime &dt2) 
{
  julnum = dt2.julnum; 
  tsec = dt2.tsec; 
  return *this;
}

HydDateTime& HydDateTime::operator=(const HydDate &dt2) 
{
  julnum = dt2.julnum; 
  return *this;
}

HydDateTime operator+ (int sec, const HydDateTime &dt) 
{
  if (sec < 0)
    return operator-(dt, abs(sec));
  else
    return HydDateTime(dt.julnum + (dt.tsec + sec) / 86400,
		    (dt.tsec + sec) % 86400);
}

HydDateTime operator+ (const HydDateTime &dt, int sec) 
{ 
  return operator+(sec, dt);
}

void operator+= (HydDateTime &dt, int sec) 
{
  if (sec < 0)
    dt -= abs(sec);
  else 
    {
      dt.julnum += ((dt.tsec + sec) / (86400)) ;
      dt.tsec = ((dt.tsec + sec) % (86400));
    }
}

HydDateTime HydDateTime::Floor(int precision_min)
{
  HydDateTime dt2(*this);
  int min=dt2.tsec/60;

  min -= min % precision_min;
  
  dt2.tsec = min*60;
  return dt2;
}

HydDateTime HydDateTime::Ceiling(int precision_min)
{
  HydDateTime dt2(*this);
  int min=dt2.tsec/60;

  if(min % precision_min)
    {
      min -= min % precision_min;
      min += precision_min;

      if(min==24*60)
	{
	  min=0;
	  dt2+=24*60*60;
	}

      dt2.tsec = min*60;
    }
  
  return dt2;
}

unsigned int operator- (const HydDateTime &dt1, const HydDateTime &dt2)
{ 
  // this func. returns the absolute difference between the two times
  // in seconds
  
  unsigned int result;
  if (dt1 >= dt2) 
    { 
      if (dt1.tsec < dt2.tsec) 
	{
	  result = dt1.tsec + 86400 - dt2.tsec;
	  result += (dt1.julnum-1 - dt2.julnum)*86400;
	}
      else 
	{
	  result = dt1.tsec - dt2.tsec;
	  result += (dt1.julnum - dt2.julnum)*86400;
	}
    }
  else 
    {
      if (dt2.tsec < dt1.tsec) 
	{
	  result = dt2.tsec + 86400 - dt1.tsec;
	  result += (dt2.julnum-1 - dt1.julnum)*86400;
	}
      else 
	{
	  result = dt2.tsec - dt1.tsec;
	  result += (dt2.julnum - dt1.julnum)*86400;
	}
    }
  return result;
} // end operator-


HydDateTime operator- (const HydDateTime &dt, int sec) 
{ 
  if (sec < 0)
    return operator+(dt, abs(sec));
  else 
    {
      if (dt.tsec < (totalSecondsType) sec) 
	{
	  julType j; totalSecondsType ts;
	  int needed = (sec - dt.tsec - 1) / (86400) + 1;
	  j = dt.julnum - needed;
	  ts = dt.tsec + (needed * 86400) - sec;
	  return HydDateTime(j, ts); 
	} 
      else 
	return HydDateTime(dt.julnum, dt.tsec - sec); 
    }
}

HydDateTime operator- (int sec, const HydDateTime &dt) 
{
  return operator-(dt, sec);
}

void operator-= (HydDateTime &dt, int sec) 
{ 
  if (sec < 0)
    dt += abs(sec);
  else 
    {
      if (dt.tsec < (totalSecondsType) sec) 
	{
	  int needed = ((sec - dt.tsec - 1) / (86400) + 1);
	  dt.julnum -= needed;
	  dt.tsec = dt.tsec + (needed * 86400) - sec;
	}
      else 
	dt.tsec -= sec; 
    }
}

int operator<(const HydDateTime &dt1, const HydDateTime &dt2)  
{ 
  return int(dt1.julnum < dt2.julnum || 
	     (dt1.julnum == dt2.julnum && dt1.tsec < dt2.tsec));
}

int operator<=(const HydDateTime &dt1, const HydDateTime &dt2) 
{ 
  return int(dt1 < dt2 || dt1 == dt2);
}


int HydDateTime::isNull() const 
{
  return (*this == NoHydDateTime);
}

ostream& operator<<(ostream &os, const HydDateTime &dt) 
{
  if (dt.isNull())
    os << "NULL";
  else {
    const HydDate& d = dt;
    const CTime& t = dt;
    os << const_cast<HydDate&>(d) << const_cast<CTime&>(t);
  }
  return os;
}

void HydDateTime::Print() 
{
  char *pt = syCh();
#ifdef MAIN
  cout << pt << endl;
#else
  Rcout << pt << std::endl;
#endif // MAIN
  delete pt;
}


void HydDateTime::getHydDateTime(yearType &yyyy, monthType &mm, dayType &dd, 
			   hourType &hh, minuteType &mn, secondType &ss) 
{
  jday2date(dd, mm, yyyy);
  totSec2time(hh, mn, ss);
}


char *HydDateTime::syCh(int form) const 
{
  char *pt = new char[100];


  if(*this == NoHydDateTime)
    snprintf(pt, 99, "null");
  else 
    {
      static yearType y = 0; static monthType m = 0; static dayType d = 0;
      static hourType h = 0; static minuteType min = 0; static secondType s = 0;
      jday2date(d, m, y);
      totSec2time(h, min, s);
      switch(form) 
	{
	case 1:
	  snprintf(pt, 99, "%02d.%02d.%d %02d:%02d", d, m, y, h, min);
	  break;
	case 2:
	  snprintf(pt, 99, "%d%02d%02d%02d%02d", y, m, d, h, min);
	  break;
	case 3:
	  snprintf(pt, 99, "%d.%02d.%02d %02d:%02d", y, m, d, h, min);
	  break;
	case 4:
	  snprintf(pt, 99, "%02d.%02d.%d",d,m,y);
	  break;
	case 5:
	  snprintf(pt, 99, "%04d%02d%02d/%02d%02d",y,m,d,h,min);
	  break;
	case 6:
	  snprintf(pt, 99, "%02d.%02d %02d:%02d",d,m,h,min);
	  break;
	case 7:
	  snprintf(pt, 99, "%02d.%02d.%d %02d:%02d:%02d", 
		  d, m, y, h, min, s);
	  break;
	default:
	  snprintf(pt, 99, "%d/%d/%d %d:%d:%d", m, d, y, h, min, s);
	  break;
	};
    }
  return (pt);
}


HydDateTime HydDateTime::StartOfDay() const 
{
  if(*this == NoHydDateTime)
    return NoHydDateTime;
  else 
    {
      static yearType y = 0; static monthType m = 0; static dayType d = 0;
      jday2date(d, m, y);
      HydDateTime dt(y, m, d);
      return dt;
    }
}

HydDateTime HydDateTime::StartOfMonth() const 
{
  if(*this == NoHydDateTime)
    return NoHydDateTime;
  else 
    {
      static yearType y = 0; static monthType m = 0; static dayType d = 0;
      jday2date(d, m, y);
      HydDateTime dt(y, m, (dayType)1);
      return dt;
    } 
}

HydDateTime HydDateTime::StartOfYear() const 
{
  if(*this == NoHydDateTime)
    return NoHydDateTime;
  else 
    {
      static yearType y = 0; static monthType m = 0; static dayType d = 0;
      jday2date(d, m, y);
      HydDateTime dt(y, (monthType)1, (dayType)1);
      return dt;
    }
}

int HydDateTime::difference_minutes(HydDateTime &dt2)
{
  int daydiff=date2jday()-dt2.date2jday();
  int hourdiff=getHour()-dt2.getHour();
  int mindiff=getMinute()-dt2.getMinute();
  int diff = daydiff*24l*60l + hourdiff*60l + mindiff;

  return diff;
}

HydDateTime  &HydDateTime::add_minutes(int num_minutes)
{  
  totalSecondsType secs = 60*((tsec/60 + num_minutes) % 1440);
  julType days = julnum + (tsec/60 + num_minutes) / 1440;

  if(num_minutes<0)
    {
      num_minutes=-num_minutes;
      secs = 60*((totalSecondsType) 
		 (((int) tsec/60 + num_minutes) % 1440));
      days=julnum + tsec/86400 - num_minutes / 1440;
    }  

  HydDateTime dt(days, secs); 
  
  *this=dt;

  return *this;
}



HydDateTime HydDateTime::EndOfMonth() const 
{
  if(*this == NoHydDateTime)
    return NoHydDateTime; 
  else 
    {    
      static yearType y = 0; static monthType m = 0; static dayType d = 0;
      jday2date(d, m, y);
      if((int)m==12)
	return EndOfYear();
      
      HydDateTime dt(y, m+1, 1);
      dt-=60;
      return dt;
    }
}

HydDateTime HydDateTime::EndOfDay() const 
{
  if(*this == NoHydDateTime)
    return NoHydDateTime;
  else 
    {
      static yearType y = 0; static monthType m = 0; static dayType d = 0;
      jday2date(d, m, y);
      HydDateTime dt(y, m, d, 23, 59);
      return dt;
    }
}

HydDateTime HydDateTime::EndOfYear() const 
{
  if(*this == NoHydDateTime)
    return NoHydDateTime;
  else 
    {
      static yearType y = 0; static monthType m = 0; static dayType d = 0;
      jday2date(d, m, y);
      HydDateTime dt(y, 12, 31, 23, 59);
      return dt;
    }
}

double HydDateTime::as_floating_point_year(void)
{
  HydDateTime start=(*this).StartOfYear();
  int sec_from_start=(*this)-start;
  int sec_in_year=86400*daysInYear();
  
  double year=getYear();
  
  year+=((double) sec_from_start)/((double) sec_in_year);
  
  return year;
}



int overlap(const HydDateTime &start1, const HydDateTime &end1,  
	    const HydDateTime &start2, const HydDateTime &end2)
{
  if(start1!=NoHydDateTime && end2!=NoHydDateTime && start1>end2)
    return 0;
  
  if(end1!=NoHydDateTime && start2!=NoHydDateTime && end1<start2)
    return 0;

  return 1;
}



// ######################################################################
// Return    : int
// Parameters: i, j - two HydDateTimes which will be compared
// Purpose   : compare routine for qsort().
// ######################################################################
int compare_datetime(const void *i, const void *j) 
{
  HydDateTime *dt1=(HydDateTime *) i;
  HydDateTime *dt2=(HydDateTime *) j;
  
  if((dt1->date2jday() > dt2->date2jday()) || 
     (dt1->date2jday() == dt2->date2jday() && 
      dt1->secofday() > dt2->secofday()))
    return -1;
  else if(dt1->date2jday() == dt2->date2jday() && 
	  dt1->secofday()  == dt2->secofday())
    return 0;
  else
    return 1;
} /* compar */

int compare_datetime_descending(const void *i, const void *j) 
{
  HydDateTime *dt1=(HydDateTime *) i;
  HydDateTime *dt2=(HydDateTime *) j;
  
  if((dt1->date2jday() > dt2->date2jday()) || 
     (dt1->date2jday() == dt2->date2jday() && 
      dt1->secofday() > dt2->secofday()))
    return 1;
  else if(dt1->date2jday() == dt2->date2jday() && 
	  dt1->secofday()  == dt2->secofday())
    return 0;
  else
    return -1;
} /* compar */

HydDateTime *find_unique_date_times(// input:
				 HydDateTime **alltimes, int *len, int num_arrays,
				 // output:
				 int *common_len)
{
  int i,j,k,totlen=0;
  for(i=0;i<num_arrays;i++)
    totlen+=len[i];

  // Make single array containing all date times:
  HydDateTime *dt1=new HydDateTime[totlen];
  k=0;
  for(i=0;i<num_arrays;i++)
    for(j=0;j<len[i];j++)
      dt1[k++]=alltimes[i][j];
  
  // sort it:
  qsort(dt1, size_t(totlen), sizeof(HydDateTime), compare_datetime);
  
  // Remove HydDateTime duplicates.

  // count number of non-duplicates:
  k=0;
  for(i=totlen-1;i>=0;i--)
    if(i==(totlen-1) || dt1[i]!=dt1[i+1])
      k++;
  
  if(k==0)
    {
      *common_len=0;
      delete [] dt1;
      return NULL;
    }

  int retlen=k;
  HydDateTime *dt2=new HydDateTime[retlen];
  k=0;
  for(i=totlen-1;i>=0;i--)
    if(i==(totlen-1) || dt1[i]!=dt1[i+1])
      dt2[k++]=dt1[i];
  
  delete [] dt1;
  
  *common_len=retlen;
  return dt2;
}





// ********************************************************
// Double linked list and file reading algorithms:
// ********************************************************


class double_linked_list
{
private:
  double_linked_list *next; // pointer to next item in a list
  double_linked_list *prev; // pointer to the previos item in a list

public: 
  // don't use these lowlevel functions unless 
  // it's absolutely neccessary!
  void setprev(double_linked_list *item);
  void setnext(double_linked_list *item);

public:
  // constructors;
  double_linked_list();     // create a period editor without surroundings
  double_linked_list(double_linked_list *head); 
  // instance that will point to the previous head of a list

  double_linked_list(double_linked_list *previtem,
		     double_linked_list *nextitem); // put the instance
                                            // in the midle of the list
  // destructor, deletes down in the link chain
  virtual ~double_linked_list();

  // Put a new element into a list with prev=prev and next=prev->next,
  // so that the new list will be prev->new->next
  void put_into(double_linked_list *prevelem);

  void removefromlist(void);        // removes this item from the list
  double_linked_list *getnext(void) {return next;} //get the next item
  double_linked_list *getprev(void) {return prev;} // get the previous item
  double_linked_list *getlast(void);       // get the last item in the list
  double_linked_list *gethead(void);       // get the first element
  // Step through the next element 'element_number' times, starting
  // from the present element , return the resulting element
  double_linked_list *get_rel_next(int elements_number);
  // Step through the next element 'element_number' times, starting
  // from the head element, return the resulting element
  double_linked_list *get_from_start(int elements_number);
  // Step through the previous element 'element_number' times, starting
  // from the present element , return the resulting element
  double_linked_list *get_rel_prev(int elements_number);
  // Step through the previous element 'element_number' times, starting
  // from the last element, return the resulting element
  double_linked_list *get_from_end(int elements_number);

  // return the number of elements from the present and to the end
  int number_of_elements(void);

  bool appears_in_list(double_linked_list *elem);

  // append one list to the end of this
  void join(double_linked_list *head2);    
  void append(double_linked_list *item); //just another name for it...

  // put a single element into the start of a list
  // update the start pointer
  friend void into_the_start(double_linked_list *new_elem,
			     double_linked_list **start);
  // put a single element into the start of a list
  // update the start and end pointers
  friend void into_the_start(double_linked_list *new_elem,
			     double_linked_list **start,
			     double_linked_list **end);
  // put a single element into the end of a list
  // update the start pointer
  friend void into_the_tail(double_linked_list *new_elem,
			    double_linked_list **start);
  // put a single element into the end of a list
  // update the start and end pointers
  friend void into_the_tail(double_linked_list *new_elem,
			    double_linked_list **start,
			    double_linked_list **end);
};

class list_2d : public double_linked_list
{
  double_2d content;
public:
  list_2d(double_2d content_, list_2d *prev);
  ~list_2d();

  double_2d get_content(void);
  void set_content(double_2d new_content);

  list_2d *suc();
  double_2d *get_array(unsigned int *len);
};

// Read a file containing two doubles per line and return
// as an array. Uses list_2d internally. Sorts the array
// according to the first double.
double_2d *get_2d_data_file(char *filename,
			    unsigned int *len, bool dosort=true);

class csv : public double_linked_list
{
  int len;
  double *val;
  char **names;

public:

  csv();
  csv(char *infile, char sep=';'); // read contents from a file (multiple rows)
  csv(csv *prev, int new_len, double *new_values); // make and set the contents of one row
  csv(csv *prev, int new_len, double *new_values, char **new_names); 
  // make and set the contents of one row
  
  ~csv();
  
  csv *prev(void);
  csv *suc(void);
  
  int get_length(void);
  double get_value(int number);
  double *get_values(void);
  // returns values not only from this element
  // but also all subsequent elements:
  double **get_all_values(int *num_elem=NULL); 
  char *get_name(int number);
  char **get_names(void);
  
  void show(ostream &out);
};

#define divfunc_h
#define lists_h

// **********************************************************************
// Differnet Aggregation Methods 
// **********************************************************************
enum METHOD 
{
  UNKNOWN_METHOD                 = -1,
  INSTANTANEOUS                  = 0,
  MAX                            = 1,
  MIN                            = 2,
  MEAN                           = 3,
  CHANGE                         = 4,
  SUM                            = 5,
  TIME_INDEPENDENT_INSTANTANEOUS = 6,
  MORNING                        = 7,
  NOON                           = 8,
  AFTERNOON                      = 9,
  MIXED_METHODS                  = 10,
  SEVERAL_YEAR_MEAN              = 11,
  SEVERAL_YEAR_MAX               = 12,
  SEVERAL_YEAR_MIN               = 13,
  SEVERAL_YEAR_SUM               = 14,
  VARIATION                      = 20,
  STANDARD_DEVIATION             = 21,
  AVERAGE_DEVIATION              = 22,
  SKEW                           = 24,
  CURTOSIS                       = 25,
  SEVERAL_YEAR_STANDARD_DEVIATION= 26,
  STANDARD_ERROR                 = 27,
  MEAN_PLUS_SDEV                 = 50,
  MEAN_MINUS_SDEV                = 51,
  SEVERAL_YEAR_MEAN_PLUS_SDEV    = 60,
  SEVERAL_YEAR_MEAN_MINUS_SDEV   = 61,
  MEAN_PLUS_SERR                 = 62,
  MEAN_MINUS_SERR                = 63,
  MEAN_PLUS_2SERR                = 64,
  MEAN_MINUS_2SERR               = 65,
  SEVERAL_YEAR_PERCENTILE_2_5    = 102,
  SEVERAL_YEAR_PERCENTILE_5      = 105,
  SEVERAL_YEAR_PERCENTILE_10     = 110,
  SEVERAL_YEAR_PERCENTILE_20     = 120,
  SEVERAL_YEAR_PERCENTILE_25     = 125,
  SEVERAL_YEAR_PERCENTILE_30     = 130,
  SEVERAL_YEAR_PERCENTILE_40     = 140,
  SEVERAL_YEAR_MEDIAN            = 150,
  SEVERAL_YEAR_PERCENTILE_60     = 160,
  SEVERAL_YEAR_PERCENTILE_70     = 170,
  SEVERAL_YEAR_PERCENTILE_75     = 175,
  SEVERAL_YEAR_PERCENTILE_80     = 180,
  SEVERAL_YEAR_PERCENTILE_90     = 190,
  SEVERAL_YEAR_PERCENTILE_95     = 195,
  SEVERAL_YEAR_PERCENTILE_97_5   = 198,
  PERCENTILE_2_5                 = 202,
  PERCENTILE_5                   = 205,
  PERCENTILE_10                  = 210,
  PERCENTILE_20                  = 220,
  PERCENTILE_25                  = 225,
  PERCENTILE_30                  = 230,
  PERCENTILE_40                  = 240,
  MEDIAN                         = 250,
  PERCENTILE_60                  = 260,
  PERCENTILE_70                  = 270,
  PERCENTILE_75                  = 275,
  PERCENTILE_80                  = 280,
  PERCENTILE_90                  = 290,
  PERCENTILE_95                  = 295,
  PERCENTILE_97_5                = 298,
  DERIVATIVE                     = 300,
  DIFFERENCE_MET                 = 301,
  ACCUMULATION_DIFFERENCE        = 302,
  CONFORM_TRAN_DAY               = 401,
  CONFORM_TRAN_WEEK              = 402,
  CONFORM_TRAN_YEAR              = 403,
  MOST_FREQUENT                  = 500
};

#define NUMBER_OF_METHODS 68

// Return true if there's enough non-missing data:
bool enough(double *data, int len, double procent_needed);

// Full summary statistics package:
double find_statistics(double *data, int len, METHOD met, 
		       double procent_needed=100.0);

// Returns the one-step auto-correlation of a sample 'x'
double get_auto_correlation(double *x, int len);  

// Functions for finding mean and standard deviation:
double find_mean(double *data, int len, 
		 double procent_needed);
double find_stdev(double *data, int len, 
		  double procent_needed, double mean=MISSING_VALUE);

// *********************************
// Linear interpolation functions:
// *********************************

double linear(double x0, double y0, double x1, double y1, double x);
double linear(HydDateTime x0, double y0, HydDateTime x1, double y1, HydDateTime x);

// Diverse functions:

int almost_equal(double v1, double v2);

// Count the number of spaces between contents
// For instant "ajksdh   skjdfh     sasasd" would return "2".
int count_spaces(char *line);

// fetches the next double and returns a pointer to 
// the end of the fetched text
char *getnextdouble(char *instr, double *nextdouble);

// fetches the next integer and returns a pointer to 
// the end of the fetched text
char *getnextint(char *instr,int *nextint); 

// fetches the next string and returns a pointer to 
// the end of the fetched text
char *getnextstr(char *instr,char *nextstr);

/* randomize using HydDateTime::now and pid */ 
int randify(bool reseed=false);

// Return random double rpecision number (strictly) 
// between 0 and 1, using "rand" in stdlib
double drand(void);

// ****************************************************
// Complex values 
// ****************************************************

class Complex
{
protected:
  double re,im;
  
public:
  
  Complex();
  Complex(double realpart);
  Complex(double realpart, double imagpart);
  Complex(double radius, double angle, bool radial // alternative, degrees 
  );
  Complex(const Complex &orig);
  
  // Setting operators
  Complex& operator= (const Complex &orig);
  Complex& operator= (const double &orig);
  Complex& operator= (const float &orig);
  Complex& operator= (const int &orig);

  // arithmetic operators
  friend Complex operator+ (const Complex &c1, const Complex &c2);
  friend Complex operator+ (const double &c1, const Complex &c2);
  friend Complex operator+ (const Complex &c1, const double &c2);
  friend Complex operator+ (const float &c1, const Complex &c2);
  friend Complex operator+ (const Complex &c1, const float &c2);
  friend Complex operator+ (const int &c1, const Complex &c2);
  friend Complex operator+ (const Complex &c1, const int &c2);

  friend Complex operator- (const Complex &c1, const Complex &c2);
  friend Complex operator- (const double &c1, const Complex &c2);
  friend Complex operator- (const Complex &c1, const double &c2);
  friend Complex operator- (const float &c1, const Complex &c2);
  friend Complex operator- (const Complex &c1, const float &c2);
  friend Complex operator- (const int &c1, const Complex &c2);
  friend Complex operator- (const Complex &c1, const int &c2);

  friend Complex operator* (const Complex &c1, const Complex &c2);
  friend Complex operator* (const double &c1, const Complex &c2);
  friend Complex operator* (const Complex &c1, const double &c2);
  friend Complex operator* (const float &c1, const Complex &c2);
  friend Complex operator* (const Complex &c1, const float &c2);
  friend Complex operator* (const int &c1, const Complex &c2);
  friend Complex operator* (const Complex &c1, const int &c2);

  friend Complex operator/ (const Complex &c1, const Complex &c2);
  friend Complex operator/ (const double &c1, const Complex &c2);
  friend Complex operator/ (const Complex &c1, const double &c2);
  friend Complex operator/ (const float &c1, const Complex &c2);
  friend Complex operator/ (const Complex &c1, const float &c2);
  friend Complex operator/ (const int &c1, const Complex &c2);
  friend Complex operator/ (const Complex &c1, const int &c2);

  friend void operator+= (Complex &c1,const Complex &c2);
  friend void operator+= (Complex &c1,const double &c2);
  friend void operator+= (Complex &c1,const float &c2);
  friend void operator+= (Complex &c1,const int &c2);

  friend void operator-= (Complex &c1,const Complex &c2);
  friend void operator-= (Complex &c1,const double &c2);
  friend void operator-= (Complex &c1,const float &c2);
  friend void operator-= (Complex &c1,const int &c2);

  friend void operator*= (Complex &c1,const Complex &c2);
  friend void operator*= (Complex &c1,const double &c2);
  friend void operator*= (Complex &c1,const float &c2);
  friend void operator*= (Complex &c1,const int &c2);

  friend void operator/= (Complex &c1,const Complex &c2);
  friend void operator/= (Complex &c1,const double &c2);
  friend void operator/= (Complex &c1,const float &c2);
  friend void operator/= (Complex &c1,const int &c2);

  friend void operator++ (Complex &c1);
  friend void operator++ (Complex &c1, int);
  friend void operator-- (Complex &c1);
  friend void operator-- (Complex &c1, int);
  
  // logical operators
  friend int operator== (const Complex &c1, const Complex &c2);
  friend int operator!= (const Complex &c1, const Complex &c2);

  double Real_value(void);
  double Imaginary_value(void);
  double Re(void);
  double Im(void);

  void set(double real, double imag);
  void set_polar(double radius, double angle);

  double angle(void);
  double radius(void);
  double abs(void);
  double abs2(void); // R^2

  Complex conjugate(void);
  Complex inverse(void); 
  Complex negative(void);

  friend double complex_abs(Complex &orig);
  friend double complex_abs2(Complex &orig);
  friend Complex complex_exp(Complex &a);
  friend Complex complex_log(Complex &a);
  friend Complex complex_sqrt(Complex &a);
  
  friend ostream& operator<<(ostream &, Complex &);
};



// ****************************************************
// Basic methods for extracting statistics 
// ****************************************************

double *mean_of_vectors(int number_of_samples, double **vectors, 
			int dimension);

double **estimated_variance(int number_of_samples,
			    double **samples /* a sample of vectors */,
			    int dimension);
double **estimated_variance(int number_of_samples,
			    double **samples /* a sample of vectors */,
			    double *mean, int dimension);


// ****************************************************
// Basic methods for sampling 
// ****************************************************

// Sample from the normal distribution:
double get_random_gauss(void);

// Sample from the multinormal distribution, using Cholesky:
double **multinormal_sample(int number_of_samples, double *expectation, 
			    double **sigma, int dimension, 
			    bool dont_keep_cholesky=false);

// ****************************************************
// Basic CDFs
// ****************************************************

double standard_normal_cdf(double x, double mean, double sigma);
double normal_cdf(double x, double mean, double sigma);
double normal_invcdf(double p, double mean, double sigma);
double chisq_deg1_cdf(double x);

// ****************************************************
// Multivariate PDFs 
// ****************************************************

double pdf_multinormal(double *x, double *expectation, double **sigma, int dimension,
		       bool logarithmic);
double pdf_multinormal(double *x, double *expectation, 
		       int dimension, double **inv_sigma, double det,
		       bool logarithmic);

double pdf_multi_t(double *x, double *expectation, double **sigma, double nu,
		   int dimension,
		   bool logarithmic);

double pdf_multi_t(double *x, double *expectation, double nu,
		   int dimension, double **inv_sigma, double det,
		   bool logarithmic);

// ****************************************************
// Linear algebra methods, taken from linalg in
// hydrasub-library.
// ****************************************************

double **Make_matrix(int rows, int columns);
Complex **Make_Complex_matrix(int rows, int columns);


// invert a given matrix;

#ifdef __cplusplus
extern "C" {
#endif
double **matrix_inverse(double **X, int dimension);
#ifdef __cplusplus
}
#endif

Complex **Complex_matrix_inverse(Complex **matrix, int dimension);


// Eigenvector decomposition:
Complex **matrix_eigenvectors(double **A, int size, Complex **eigen_values,
			      Complex ***V_inv=NULL);
Complex *Complex_eigenvalues(double **A, int size);
double *double_eigenvalues(double **A, int size);

// find the determinant of a given matrix using eigenvector analysis;
double matrix_determinant(double **matrix_, int dimension);
Complex log_matrix_determinant(double **matrix_, int dimension);

// Cholesky decomposition:
double **get_cholesky(double **matrix_, int dimension);






// ****************************************************
// Definitions done
// ****************************************************



// Complex class:

Complex::Complex()
{
  re=im=0.0;
}

Complex::Complex(double realpart)
{
  re=realpart;
  im=0.0;
}

Complex::Complex(double realpart, double imagpart)
{
  re=realpart;
  im=imagpart;
}

Complex::Complex(double radius, double angle, bool radial // alternative, degrees 
		   )
{
  double angle2=radial ? angle : 2.0*M_PI*angle/360.0;
  re=radius*cos(angle2);
  im=radius*sin(angle2);
}

Complex::Complex(const Complex &orig)
{
  re=orig.re;
  im=orig.im;
}
  
// Setting operators
Complex& Complex::operator= (const Complex &orig)
{
  re=orig.re;
  im=orig.im;
  return *this;
}

Complex& Complex::operator= (const double &orig)
{
  re=orig;
  im=0.0;
  return *this;
}

Complex& Complex::operator= (const float &orig)
{
  re=(double) orig;
  im=0.0;
  return *this;
}

Complex&Complex:: operator= (const int &orig)
{
  re=(double) orig;
  im=0.0;
  return *this;
}

Complex operator+ (const Complex &c1, const Complex &c2)
{
  double re=c1.re+c2.re, im=c1.im+c2.im;
  return Complex(re, im);
}

Complex operator+ (const Complex &c1, const double &c2)
{
  double re=c1.re+c2, im=c1.im;
  return Complex(re, im);
}

Complex operator+ (const double &c1, const Complex &c2)
{
  double re=c1+c2.re, im=c2.im;
  return Complex(re, im);
}

Complex operator+ (const Complex &c1, const float &c2)
{
  double re=c1.re+double(c2), im=c1.im;
  return Complex(re, im);
}

Complex operator+ (const float &c1, const Complex &c2)
{
  double re=double(c1)+c2.re, im=c2.im;
  return Complex(re, im);
}

Complex operator+ (const Complex &c1, const int &c2)
{
  double re=c1.re+double(c2), im=c1.im;
  return Complex(re, im);
}

Complex operator+ (const int &c1, const Complex &c2)
{
  double re=double(c1)+c2.re, im=c2.im;
  return Complex(re, im);
}

Complex operator- (const Complex &c1, const Complex &c2)
{
  double re=c1.re-c2.re, im=c1.im-c2.im;
  return Complex(re, im);
}

Complex operator- (const Complex &c1, const double &c2)
{
  double re=c1.re-c2, im=c1.im;
  return Complex(re, im);
}

Complex operator- (const double &c1, const Complex &c2)
{
  double re=c1-c2.re, im=-c2.im;
  return Complex(re, im);
}

Complex operator- (const Complex &c1, const float &c2)
{
  double re=c1.re-double(c2), im=c1.im;
  return Complex(re, im);
}

Complex operator- (const float &c1, const Complex &c2)
{
  double re=double(c1)-c2.re, im=-c2.im;
  return Complex(re, im);
}

Complex operator- (const Complex &c1, const int &c2)
{
  double re=c1.re-double(c2), im=c1.im;
  return Complex(re, im);
}

Complex operator- (const int &c1, const Complex &c2)
{
  double re=double(c1)-c2.re, im=-c2.im;
  return Complex(re, im);
}

Complex operator* (const Complex &c1, const Complex &c2)
{
  double re=c1.re*c2.re-c1.im*c2.im, im=c1.im*c2.re+c1.re*c2.im;
  return Complex(re, im);
}

Complex operator* (const Complex &c1, const double &c2)
{
  double re=c1.re*c2, im=c1.im*c2;
  return Complex(re, im);
}

Complex operator* (const double &c1, const Complex &c2)
{
  double re=c1*c2.re, im=c1*c2.im;
  return Complex(re, im);
}

Complex operator* (const Complex &c1, const float &c2)
{
  double re=c1.re*double(c2), im=c1.im*double(c2);
  return Complex(re, im);
}

Complex operator* (const float &c1, const Complex &c2)
{
  double re=double(c1)*c2.re, im=double(c1)*c2.im;
  return Complex(re, im);
}

Complex operator* (const Complex &c1, const int &c2)
{
  double re=c1.re*double(c2), im=c1.im*double(c2);
  return Complex(re, im);
}

Complex operator* (const int &c1, const Complex &c2)
{
  double re=double(c1)*c2.re, im=double(c1)*c2.im;
  return Complex(re, im);
}

Complex operator/ (const Complex &c1, const Complex &c2)
{
  double re=(c1.re*c2.re+c1.im*c2.im)/(c2.re*c2.re+c2.im*c2.im);
  double im=(c1.im*c2.re-c1.re*c2.im)/(c2.re*c2.re+c2.im*c2.im);
  return Complex(re, im);
}

Complex operator/ (const Complex &c1, const double &c2)
{
  double re=c1.re/c2, im=c1.im/c2;
  return Complex(re, im);
}

Complex operator/ (const double &c1, const Complex &c2)
{
  double re=c1*c2.re/(c2.re*c2.re+c2.im*c2.im);
  double im=c1*c2.im/(c2.re*c2.re+c2.im*c2.im);
  return Complex(re, im);
}

Complex operator/ (const Complex &c1, const float &c2)
{
  double re=c1.re/double(c2), im=c1.im/double(c2);
  return Complex(re, im);
}

Complex operator/ (const float &c1, const Complex &c2)
{
  double re=double(c1)*c2.re/(c2.re*c2.re+c2.im*c2.im);
  double im=double(c1)*c2.im/(c2.re*c2.re+c2.im*c2.im);
  return Complex(re, im);
}

Complex operator/ (const Complex &c1, const int &c2)
{
  double re=c1.re/double(c2), im=c1.im/double(c2);
  return Complex(re, im);
}

Complex operator/ (const int &c1, const Complex &c2)
{
  double re=double(c1)*c2.re/(c2.re*c2.re+c2.im*c2.im);
  double im=double(c1)*c2.im/(c2.re*c2.re+c2.im*c2.im);
  return Complex(re, im);
}

void operator+= (Complex &c1,const Complex &c2)
{
  c1=c1+c2;
}

void operator+= (Complex &c1,const double &c2)
{
  c1=c1+c2;
}

void operator+= (Complex &c1,const float &c2)
{
  c1=c1+c2;
}

void operator+= (Complex &c1,const int &c2)
{
  c1=c1+c2;
}

void operator-= (Complex &c1,const Complex &c2)
{
  c1=c1-c2;
}

void operator-= (Complex &c1,const double &c2)
{
  c1=c1-c2;
}

void operator-= (Complex &c1,const float &c2)
{
  c1=c1-c2;
}

void operator-= (Complex &c1,const int &c2)
{
  c1=c1-c2;
}

void operator*= (Complex &c1,const Complex &c2)
{
  c1=c1*c2;
}

void operator*= (Complex &c1,const double &c2)
{
  c1=c1*c2;
}

void operator*= (Complex &c1,const float &c2)
{
  c1=c1*c2;
}

void operator*= (Complex &c1,const int &c2)
{
  c1=c1*c2;
}

void operator/= (Complex &c1,const Complex &c2)
{
  c1=c1/c2;
}

void operator/= (Complex &c1,const double &c2)
{
  c1=c1/c2;
}

void operator/= (Complex &c1,const float &c2)
{
  c1=c1/c2;
}

void operator/= (Complex &c1,const int &c2)
{
  c1=c1/c2;
}

void operator++ (Complex &c1)
{
  c1=c1+1;
}
  
void operator++ (Complex &c1, int)
{
  c1=c1+1;
}

void operator-- (Complex &c1)
{
  c1=c1-1;
}

void operator-- (Complex &c1, int)
{
  c1=c1-1;
}

// logical operators
int operator== (const Complex &c1, const Complex &c2)
{
  return int(c1.re==c2.re && c1.im==c2.im);
}

int operator!= (const Complex &c1, const Complex &c2)
{
  return int(c1.re!=c2.re || c1.im!=c2.im);
}

double Complex::Real_value(void)
{
  return re;
}

double Complex::Imaginary_value(void)
{
  return im;
}

double Complex::Re(void)
{
  return re;
}

double Complex::Im(void)
{
  return im;
}

void Complex::set(double real, double imag)
{
  re=real;
  im=imag;
}

void Complex::set_polar(double radius, double angle)
{
  re=radius*cos(angle);
  im=radius*sin(angle);
}


double Complex::angle(void)
{
  return atan(im/re);
}

double Complex::radius(void)
{
  return sqrt(re*re+im*im);
}

double Complex::abs(void)
{
  return sqrt(re*re+im*im);
}

double Complex::abs2(void) // R^2
{
  return re*re+im*im;
}

Complex Complex::conjugate(void)
{
  Complex ret(re,-im);
  return ret;
}

Complex Complex::inverse(void)
{
  Complex ret=1.0/(*this);
  return ret;
}

Complex Complex::negative(void)
{
  Complex ret(-re,-im);
  return ret;
}


Complex complex_exp(Complex &a)
{
  Complex ret(exp(a.Re()),a.Im(),true);
  return ret;
}


Complex complex_log(Complex &a)
{
  Complex ret(0.5*log(a.Re()*a.Re()+a.Im()*a.Im()), atan(a.Im()/a.Re()));
  return ret;
}

Complex complex_sqrt(Complex &a)
{
  Complex ret(sqrt(a.radius()),a.angle()/2.0,true);
  return ret;
}


ostream& operator<<(ostream &out,  Complex &c)
{
  char line[100];
  double re=c.Re();
  double im=c.Im();

  snprintf(line, 99, "%8.4g+%8.4gi", re,im);
  out << line;

  return out;
}


// Complex class done





// Count the number of spaces between contents
// For instant "ajksdh   skjdfh     sasasd" would return "2".
int count_spaces(char *line)
{
  int count=0;
  int has_been_content=0, has_been_spaces=0;

  for(unsigned int i=0;i<strlen(line);i++)
    {
      if(line[i]!=' ' && line[i]!='\t')
	{
	  has_been_content=1;
	  if(has_been_spaces)
	    count++;
	  has_been_spaces=0;
	}
      else if(has_been_content)
	has_been_spaces=1;
    }

  return count;
}


// fetches the next double and returns a pointer to 
// the end of the fetched text
char *getnextdouble(char *instr, double *nextdouble)
{
  char *ptr=instr;
  
  *nextdouble=0.0;

  while(!isdigit(*ptr) && *ptr!='+' && *ptr!='-' &&
	strncmp(ptr,"---",3) && strncasecmp(ptr, "missing", 7) &&
	strncasecmp(ptr, "na", 2))
    {
      if(!(*ptr))
	return NULL;
      
      ptr++;
    }

  
  if(!strncmp(ptr,"---",3) || !strncasecmp(ptr, "missing", 7) || 
     !strncasecmp(ptr, "na", 2))
    *nextdouble=MISSING_VALUE;
  else if(!sscanf(ptr,"%lf",nextdouble))
    return NULL;
  
  while(*ptr && *ptr!=' ' && *ptr!='\t' && *ptr!=',' && *ptr!=';' && *ptr!=':')
    ptr++;
  
  return ptr;
}


// fetches the next integer and returns a pointer to 
// the end of the fetched text
char *getnextint(char *instr,int *nextint)
{
  char *ptr=instr;
  
  *nextint=0;

  while(!isdigit(*ptr) && *ptr!='-')
    {
      if(!(*ptr))
	return NULL;
      
      ptr++;
    }
  
  if(!sscanf(ptr,"%d",nextint))
    return NULL;
  
  while(*ptr && isdigit(*ptr))
    ptr++;

  return ptr;
}


// fetches the next string and returns a pointer to 
// the end of the fetched text
char *getnextstr(char *instr,char *nextstr)
{
  int hasdelimiter=0;
  char *ptr=instr;
  
  strcpy(nextstr,"");

  while(*ptr==' ' || *ptr==',' || *ptr=='\t')
    {
      if(!(*ptr))
	return NULL;
      
      if(*ptr=='\"')
	{
	  hasdelimiter=1;
	  ptr++;
	  break;
	}
      
      ptr++;
    }

  
  char *nextptr=nextstr;

  while(*ptr && (hasdelimiter || (*ptr!=' ' && *ptr!=',' && *ptr!='\t')))
    {
      if(hasdelimiter && *ptr=='\"')
	{
	  ptr++;
	  break;
	}
      
      *nextptr=*ptr;
      ptr++;
      nextptr++;
    }
  *nextptr='\0';

  return ptr;
}


int randify(bool reseed)
{
  static bool visited=false;
  HydDateTime now,then(1970,1,1); 
  now.now();
  
  long int seed=(now-then);
  
  if(!visited || reseed)
    {
#ifdef MAIN
      seed+=getpid();
      srand(seed);
#else
      //set_seed(seed);
#endif // MAIN
    }
      
  visited=true;
  
  return seed;
}


double drand(void)
{
  double ret=0.5;

  do
    {
#ifdef MAIN
      int r=rand();
      int rmax=RAND_MAX;
      
      if(rmax<100000)
	ret=double(r*rmax+rand())/double(rmax*rmax);
      else
	ret=double(r)/double(rmax);
#else
      ret=unif_rand();
#endif // MAIN      
    } while(ret==0.0 || ret==1.0);
  
  return ret;
}

// *******************************************************
//                   DOUBLE_LINKED_LIST
// Just as the name says, a double linked list, used
// as a parent class when implementing various classes
// that needs to be linked lists.
// *******************************************************



// constructors;

// constructor
// create a period editor without surroundings
double_linked_list::double_linked_list()
{
  next=prev=NULL;
} // double_linked_list::double_linked_list


// constructor
// object that will point to the previous head of a list
double_linked_list::double_linked_list(double_linked_list *head) 
{
  next=head;
  if(head)
    head->prev=this;
  prev=NULL;
} // double_linked_list::double_linked_list


// constructor
// put the object in the middle of the list
double_linked_list::double_linked_list(double_linked_list *previtem,
				       double_linked_list *nextitem)
{
  next=nextitem;
  prev=previtem;
  
  if(prev)
    prev->next=this;

  if(next)
    next->prev=this;
} // double_linked_list::double_linked_list


// destructor
// delete recursively down a chain...
double_linked_list::~double_linked_list()
{
  if(next)
    delete next;
} // double_linked_list::~double_linked_list


  // Put a new element into a list with prev=prevelem and next=prevelem->next,
  // so that the new lsit will be prevelem->new->next
void double_linked_list::put_into(double_linked_list *prevelem)
{
  prev=prevelem;
  if(prev)
    next=prevelem->next;
  else
    next=NULL;
  if(next)
    next->prev=this;
  if(prev)
    prev->next=this;
}


//                 REMOVEFROMLIST
// removes an object from a linked list...
void double_linked_list::removefromlist(void)
{
  if(prev)
    prev->next=next;
  if(next)
    next->prev=prev;
  prev=NULL;
  next=NULL;
} // double_linked_list::removefromlist


//                 GETLAST
// Fetches the last object in a list
double_linked_list *double_linked_list::getlast(void)
{
  double_linked_list *pt;
  for(pt=this;pt->getnext();pt=pt->getnext()); // get the last item
  return pt;
} // double_linked_list::getlast


//                 GETHEAD
// Fetches the first object in a list
double_linked_list *double_linked_list::gethead(void)
{
  double_linked_list *pt;
  for(pt=this;pt->getprev();pt=pt->getprev()); // get the last item
  return pt;
} // double_linked_list::getlast


//             GET_REL_NEXT
// Gets the element 'num' places next to this element
double_linked_list * double_linked_list::get_rel_next(int num)
{
  int i;
  double_linked_list *pt=this;
  for(i=1;i<=num && pt;i++)
    pt=pt->getnext();
  return pt;
}

//             GET_FROM_START
// Gets the element 'num' places next to the head of the list
double_linked_list *double_linked_list::get_from_start(int num)
{
  int i;
  double_linked_list *pt=this->gethead();
  for(i=1;i<=num && pt;i++)
    pt=pt->getnext();
  return pt;
}


//             GET_REL_PREV
// Gets the element 'num' places previous to this element
double_linked_list *double_linked_list::get_rel_prev(int num)
{
  int i;
  double_linked_list *pt=this;
  for(i=1;i<=num && pt;i++)
    pt=pt->getprev();
  return pt;
}

//             GET_FROM_END
// Gets the element 'num' places previous to the end of the list
double_linked_list *double_linked_list::get_from_end(int num)
{
  int i;
  double_linked_list *pt=this->getlast();
  for(i=1;i<=num && pt;i++)
    pt=pt->getprev();
  return pt;
}

bool double_linked_list::appears_in_list(double_linked_list *elem)
{
  double_linked_list *pt=this;

  do
    {
      if(pt==elem)
	return true;
    } while((pt=pt->getnext()) != NULL);
  
  return false;
}

//                  JOIN
// append this list with the head of another...
void double_linked_list::join(double_linked_list *head2)
{
  double_linked_list  *pt=getlast();

  pt->next=head2;
  head2->prev=pt;
} // double_linked_list::join

// uses 'join'
void double_linked_list::append(double_linked_list *item)
{
  join(item);
}

// counts the number of elements following (and including) the present
int double_linked_list::number_of_elements(void)
{
  int i=0;
  for(double_linked_list *pt=this; pt; pt=pt->getnext())
    i++;
  return i;
}

void double_linked_list::setprev(double_linked_list *item)
{
  prev=item;
  if(item)
    item->next=this;
}

void double_linked_list::setnext(double_linked_list *item)
{
  next=item;
  if(item)
    item->prev=this;
}

// Functions used for adding elements into a double_linked_list

void into_the_start(double_linked_list *new_elem,
		    double_linked_list **start)
{
  if(!*start)
    *start=new_elem;
  else
    {
      new_elem->setnext(*start);
      *start=new_elem;
    }
}

// put a single element into the tail of a list
// update the start pointer
void into_the_tail(double_linked_list *new_elem,
		   double_linked_list **start)
{
  if(!*start)
    *start=new_elem;
  else
    {
      double_linked_list *end=(*start)->getlast();
      new_elem->setprev(end);
    }
}

// put a single element into the start of a list
// update the start and end pointers
void into_the_start(double_linked_list *new_elem,
		    double_linked_list **start,
		    double_linked_list **end)
{
  if(!*start)
    *start=*end=new_elem;
  else
    {
      new_elem->setnext(*start);
      *start=new_elem;
    }
}

// put a single element into the end of a list
// update the start and end pointers
void into_the_tail(double_linked_list *new_elem,
		   double_linked_list **start,
		   double_linked_list **end)
{
  if(!*start)
    *start=*end=new_elem;
  else
    {
      new_elem->setprev(*end);
      *end=new_elem;
    }
}


// ######################################################################
// Return    : int
// Parameters: i, j - two values which will be compared
// Purpose   : compare routine for double_2d sorting using qsort() 
// ######################################################################
int compare_double_2d_x(const void *i, const void *j) 
{
  double v =  ((double_2d *)i)->x - ((double_2d *)j)->x;
  if(v < 0)
    return -1;
  else if(v == 0)
    return 0;
  else
    return 1;
} 


// ######################################################################
// Return    : int
// Parameters: i, j - two values which will be compared
// Purpose   : compare routine for double_2d sorting using qsort() 
// ######################################################################
int compare_double_2d_y(const void *i, const void *j) 
{
  double v =  ((double_2d *)i)->y - ((double_2d *)j)->y;
  if(v < 0)
    return -1;
  else if(v == 0)
    return 0;
  else
    return 1;
} 

list_2d::list_2d(double_2d content_, list_2d *prev) :
  double_linked_list((double_linked_list *) prev, NULL)
{
  content=content_;
}

list_2d::~list_2d()
{
  list_2d *next=suc();

  this->removefromlist();

  delete next;
}

double_2d list_2d::get_content(void)
{
  return content;
}

void list_2d::set_content(double_2d new_content)
{
  content=new_content;
}

list_2d *list_2d::suc()
{
  return (list_2d *) getnext(); 
}

double_2d *list_2d::get_array(unsigned int *len)
{
  *len=number_of_elements();

  unsigned int i=0;
  double_2d *ret=new double_2d[*len];
  for(list_2d *ptr=this;ptr;ptr=ptr->suc(),i++)
    ret[i]=ptr->get_content();

  return ret;
}





double_2d *get_2d_data_file(char *filename, unsigned int *len, bool dosort)
{
  list_2d *head=NULL, *tail=NULL;
  ifstream in;
  char line[10000];

  in.open(filename, ios::in);
  if(in.fail())
    {
#ifdef MAIN
      cerr << "Couldn't open file \"" << filename << "\"" << endl;
      exit(0);
#else
      Rcout << "Couldn't open file \"" << filename << "\"" << std::endl;
#endif // MAIN
    }
  
  in.getline(line,9999);
  while(!in.eof())
    {
      double x,y;
      double_2d newcontent;

      if(sscanf(line,"%lf %lf", &x, &y)==2)
	{
	  newcontent.x=x;
	  newcontent.y=y;

	  tail=new list_2d(newcontent,tail);
	  if(!head)
	    head=tail;
	}

      in.getline(line,9999);
    }

  double_2d *ret=head->get_array(len);
  if(dosort)
    qsort(ret, (size_t) *len, sizeof(double_2d), compare_double_2d_x);
  
  delete head;

  return ret;
}






csv::csv() : double_linked_list()
{
  len=0;
  val=NULL;
  names=NULL;
}

csv::csv(char *infile, char sep) : double_linked_list()
{
  ifstream in;
  char line[10000];
  csv *head=this, *tail=NULL;

  in.open(infile, ios::in);
  if(in.fail())
    {
#ifdef MAIN
      cerr << "Could not open the file " << infile << "!" << endl;
#endif // MAIN
      return;
    }
  
  in.getline(line,9999);
  len=1;
  int i,slen=strlen(line);
  for(i=0;i<slen;i++)
    if(line[i]==sep)
      len++;
  names=new char*[len];
  for(i=0;i<len;i++)
    names[i]=new char[100];
  int j=0,k=0;
  for(i=0;i<slen;i++)
    {
      if(line[i]!=sep)
	{
	  names[j][k]=line[i];
	  k++;
	}
      
      if(line[i]==sep || i==(slen-1))
	{
	  names[j][k]='\0';
	  j++;
	  k=0;
	}
    }
  
  int lineno=1;
  in.getline(line,9999);
  while(!in.eof())
    {
      char *ptr=line;
      double *vals=new double[len];

      lineno++;

      for(i=0;i<len;i++)
	{
	  vals[i]=MISSING_VALUE;
	  if(ptr[1]==sep || ptr[1]=='\0' || ptr[1]=='\n' || ptr[1]==(char)13)
	    ptr++;
	  else
	    {
	      ptr=getnextdouble(ptr,vals+i);
	      if(i<(len-1) && ptr==NULL)
		{
#ifdef MAIN
		  cerr << i << " " << len << " " << endl; // (long int) ptr << endl;
		  cerr << "Error in line " << lineno << 
		    "! Too few values!" << endl;
		  cerr << line << endl;
#else // R context
		  Rcout << i << " " << len << " " << std::endl; // (long int) ptr << endl;
		  Rcout << "Error in line " << lineno << 
		    "! Too few values!" << std::endl;
		  Rcout << line << std::endl;
#endif // MAIN
		  if(head)
		    delete head;
		  delete [] vals;
		  return;
		}
	    }
	}
	 
      if(lineno==2)
	{
	  val=new double[len];
	  for(i=0;i<len;i++)
	    val[i]=vals[i];
	  tail=head;
	}
      else
	tail=new csv(tail, len, vals);

      in.getline(line,9999);
    }
  in.close();

  return;
}

csv::csv(csv *prev, int new_len, double *new_values) : 
  double_linked_list((double_linked_list *) prev, NULL) 
{
  len=new_len;
  val=new double[len];
  for(int i=0;i<len;i++)
    val[i]=new_values[i];
  names=NULL;
}


csv::csv(csv *prev, int new_len, double *new_values, char **new_names) : 
  double_linked_list((double_linked_list *) prev, NULL) 
{
  len=new_len;
  val=new double[len];
  names=new char*[len];
  for(int i=0;i<len;i++)
    {
      names[i]=new char[100];
      val[i]=new_values[i];
      strcpy(names[i],new_names[i]);
    }
}

csv::~csv()
{
  if(len>0)
    {
      doubledelete(names,len);
      if(val)
	delete [] val;
    }
}

csv *csv::prev(void)
{
  return (csv *) getprev();
}

csv *csv::suc(void)
{
  return (csv *) getnext();
}


int csv::get_length(void)
{
  return len;
}

double csv::get_value(int number)
{
  if(number<0 || number>len || val==NULL)
    {
#ifdef MAIN
      cerr << "Programming error - number out of range in csv::get_value" << endl;
#endif // MAIN
      return MISSING_VALUE;
    }
  return val[number];
}

double *csv::get_values(void)
{
  return val;
}

// returns values not only from this element
// but also all subsequent elements:
double **csv::get_all_values(int *num_elem)
{
  int i=0,j,num=number_of_elements();
  double **ret=new double*[num];
  
  for(csv *ptr=this;ptr;ptr=ptr->suc())
    {
      ret[i]=new double[len];
      for(j=0;j<len;j++)
	ret[i][j]=ptr->val[j];
      i++;
    }
  
  if(num_elem)
    *num_elem=num;

  return ret;
}

char *csv::get_name(int number)
{
  if(number<0 || number>len || names==NULL)
    {
#ifdef MAIN
      cerr << "Programming error - number out of range in csv::get_name" << endl;
#endif // MAIN
      return NULL;
    }
  return names[number];
}

char **csv::get_names(void)
{
  return names;
}

void csv::show(ostream &out)
{
  int i;
  for(i=0;i<len;i++)
    {
      out << names[i];
      if(i<(len-1))
	out << ";";
    }
  out << std::endl;

  for(csv *ptr=this;ptr;ptr=ptr->suc())
    {
      for(i=0;i<len;i++)
	{
	  if(ptr->val[i]!=MISSING_VALUE)
	    out << ptr->val[i];
	  if(i<(len-1))
	    out << ";";
	}
      out << std::endl;
    }
}












bool enough(double *data, int len, double procent_needed)
{
  int miss_found=0;

  for(int i=0;i<len;i++)
    if(data[i]==MISSING_VALUE)
      miss_found++;
  
  if(double(len-miss_found)*100.0 < procent_needed)
    return false;
  else
    return true;
}

// FIND_MEAN
// Finds the mean value of a timeserie
// Parameters: data - the timeserie
//             len  - the length of the serie
//             procent_needed - the procent of the data needed tp
//                              extract the mean value
double find_mean(double *data, int len, double procent_needed)
{
  double meanval=0.0;
  int num=0;
  int requiredlen=(int) (((double) len)*procent_needed/100.0);

  // loop through the data...
  for(int i=0;i<len;i++)
    if(data[i]!=MISSING_VALUE)
      {
	meanval+=data[i]; // increment the mean value
	num++; // increment the number of values used
      }

  if(num<requiredlen) // if th enumber of values found was less than needed...
    return MISSING_VALUE;

  return meanval/num;
}

// ######################################################################
// Return    : int
// Parameters: i, j - two values which will be compared
// Purpose   : compare routine for qsort() in the 'percentile' method
// ######################################################################
int compare_double(const void *i, const void *j) 
{
  double v =  *(double *)i - *(double *)j;
  if(v < 0)
    return -1;
  else if(v == 0)
    return 0;
  else
    return 1;
} /* compar */

// FIND_STDEV
// Find the standard deviation of a time serie
// Parameters: data - the serie
//             len  - the length of the serie
//             procent_needed - the procent of the data needed tp
//                              extract the mean value
//             mean - the mean value of the time serie 
//                    (can be set to MISSING_VALUE)
double find_stdev(double *data, int len, 
		  double procent_needed, double mean)
{
  double var=0.0;
  int num=0;
  int requiredlen=(int) (((double) len)*procent_needed/100.0);

  if(mean==MISSING_VALUE)
    mean=find_mean(data, len, procent_needed);

  if(mean==MISSING_VALUE)
    return MISSING_VALUE;

  // loop through the data...
  for(int i=0;i<len;i++)
    if(data[i]!=MISSING_VALUE)
      {
	// increment the variance
	var += (data[i]-mean)*(data[i]-mean); 
	
	num++; // increment the number of values used
      }

  if(num<requiredlen) // if the number of values found was less than needed...
    return MISSING_VALUE;

  if(var==MISSING_VALUE)
    return var;

  return sqrt(var/double(num-1));
}

double find_statistics(double *data, int len, METHOD met, 
		       double procent_needed)
{
  if(!enough(data, len, procent_needed))
    return MISSING_VALUE;

  if(met==MEAN || met==SEVERAL_YEAR_MEAN)
    return find_mean(data, len, procent_needed); // in order to speed things up
  else if(met==STANDARD_DEVIATION || met==SEVERAL_YEAR_STANDARD_DEVIATION)
    return find_stdev(data, len, procent_needed);
  else if(met==STANDARD_ERROR)
    return find_stdev(data, len, procent_needed)/sqrt(double(len));
  else if(met==MEAN_PLUS_SDEV || met==SEVERAL_YEAR_MEAN_PLUS_SDEV || 
	  met==MEAN_PLUS_SERR || met==MEAN_PLUS_2SERR)
    {
      double mean=find_mean(data, len, procent_needed);
      double sdev=find_stdev(data, len, procent_needed, mean);

      if(met==MEAN_PLUS_SERR)
	sdev/=sqrt(double(len));
      else if(met==MEAN_PLUS_2SERR)
	sdev*=2.0/sqrt(double(len));

      return mean+sdev;
    }
  else if(met==MEAN_MINUS_SDEV || met==SEVERAL_YEAR_MEAN_MINUS_SDEV || 
	  met==MEAN_MINUS_SERR || met==MEAN_MINUS_2SERR)
    {
      double mean=find_mean(data, len, procent_needed);
      double sdev=find_stdev(data, len, procent_needed, mean);

      if(met==MEAN_MINUS_SERR)
	sdev/=sqrt(double(len));
      else if(met==MEAN_MINUS_2SERR)
	sdev*=2.0/sqrt(double(len));

      return mean-sdev;
    }

  double val=MISSING_VALUE;
  int i;
      
  if((met>=SEVERAL_YEAR_PERCENTILE_2_5 && met<=SEVERAL_YEAR_PERCENTILE_97_5) ||
     (met>=PERCENTILE_2_5 && met<=PERCENTILE_97_5))
    {
      double *databuffer=new double[len];
      double perc=50.0;

      switch(met)
	{
	case SEVERAL_YEAR_PERCENTILE_2_5:
	case PERCENTILE_2_5:
	  perc=2.5;
	  break;
	case SEVERAL_YEAR_PERCENTILE_5:
	case PERCENTILE_5:
	  perc=5.0;
	  break;
	case SEVERAL_YEAR_PERCENTILE_10:
	case PERCENTILE_10:
	  perc=10.0;
	  break;
	case SEVERAL_YEAR_PERCENTILE_20:
	case PERCENTILE_20:
	  perc=20.0;
	  break;
	case SEVERAL_YEAR_PERCENTILE_25:
	case PERCENTILE_25:
	  perc=25.0;
	  break;
	case SEVERAL_YEAR_PERCENTILE_30:
	case PERCENTILE_30:
	  perc=30.0;
	  break;
	case SEVERAL_YEAR_PERCENTILE_40:
	case PERCENTILE_40:
	  perc=40.0;
	  break;
	case SEVERAL_YEAR_MEDIAN:
	case MEDIAN:
	  perc=50.0;
	  break;
	case SEVERAL_YEAR_PERCENTILE_60:
	case PERCENTILE_60:
	  perc=60.0;
	  break;
	case SEVERAL_YEAR_PERCENTILE_70:
	case PERCENTILE_70:
	  perc=70.0;
	  break;
	case SEVERAL_YEAR_PERCENTILE_75:
	case PERCENTILE_75:
	  perc=75.0;
	  break;
	case SEVERAL_YEAR_PERCENTILE_80:
	case PERCENTILE_80:
	  perc=80.0;
	  break;
	case SEVERAL_YEAR_PERCENTILE_90:
	case PERCENTILE_90:
	  perc=90.0;
	  break;
	case SEVERAL_YEAR_PERCENTILE_95:
	case PERCENTILE_95:
	  perc=95.0;
	  break;
	case SEVERAL_YEAR_PERCENTILE_97_5:
	case PERCENTILE_97_5:
	  perc=97.5;
	  break;
	default:
	  break;
	}
      
      int newlen=0;
      for(i=0;i<len;i++)
	if(data[i]!=MISSING_VALUE)
	  {
	    databuffer[newlen]=data[i];
	    newlen++;
	  }

      if(newlen==0)
	{
	  delete [] databuffer;
	  return MISSING_VALUE;
	}
      
      int p = int((newlen-1) * perc/100.0), q = int((newlen-1) * perc/100.0 + 
						    0.99999);

      qsort(databuffer, newlen, sizeof(double), compare_double);

      if(p == q)		// Match a perfect index
	val = databuffer[p];
      else 
	/*
	 * Lies between to data points. Percentil is calculated
	 * with linear interpolation.
	 */
	val = linear(p, databuffer[p], q, databuffer[q], 
		     double(newlen-1) * perc/100.0);
      
      delete [] databuffer;
    }
  else if(met==MAX || met==SEVERAL_YEAR_MAX)
    {
      val=data[0];
      for(i=1;i<len;i++)
	if(val==MISSING_VALUE || (data[i]>val && data[i]!=MISSING_VALUE))
	  val=data[i];
    }
  else if(met==MIN || met==SEVERAL_YEAR_MIN)
    {
      val=data[0];
      for(i=1;i<len;i++)
	if(val==MISSING_VALUE || (data[i]<val && data[i]!=MISSING_VALUE))
	  val=data[i];
    }
  else if(met==SUM || met==SEVERAL_YEAR_SUM)
    {
      val=0.0;
      for(i=0;i<len;i++)
	if(data[i]!=MISSING_VALUE)
	  val += data[i];
    }
  else if(met==VARIATION)
    {
      double mean=find_mean(data, len, procent_needed);
      int found=0;

      val=0.0;
      for(i=0;i<len;i++)
	if(data[i]!=MISSING_VALUE)
	  {
	    val += (data[i]-mean)*(data[i]-mean);
	    found++;
	  }

      val /= double(found-1);
    }
  else if(met==AVERAGE_DEVIATION)
    {
      int found=0;

      val = 0.0;
      for(i=1;i<len;i++)
	if(data[i]!=MISSING_VALUE && data[i-1]!=MISSING_VALUE)
	  {
	    // cumbersome method, but in case of missing data, we'll need it
	    val += (data[i]-data[i-1]);
	    found++;
	  }
      if(data[0]!=MISSING_VALUE && data[len-1]!=MISSING_VALUE)
	{
	  val += (data[0]-data[len-1]);
	  found++;
	}

      val /= double(found);
    }
  else if(met==SKEW)
    {
      double mean=find_mean(data, len, procent_needed);
      double sdev=find_stdev(data, len, procent_needed, mean);
      int found=0;

      val=0.0;
      for(i=0;i<len;i++)
	if(data[i]!=MISSING_VALUE)
	  {
	    val += (data[i]-mean)*(data[i]-mean)*(data[i]-mean);
	    found++;
	  }

      val *= double(found)/double(found-1)/double(found-2);
      val /= (sdev*sdev*sdev);
    }
  else if(met==CURTOSIS)
    {
      double mean=find_mean(data, len, procent_needed);
      double sdev=find_stdev(data, len, procent_needed, mean);
      int found=0;

      val=0.0;
      for(i=0;i<len;i++)
	if(data[i]!=MISSING_VALUE)
	  {
	    val += (data[i]-mean)*(data[i]-mean)*(data[i]-mean)*(data[i]-mean);
	    found++;
	  }

      val *= double(found)*double(found+1)/double(found-1)/double(found-1)/
	double(found-1);
      val /= (sdev*sdev*sdev*sdev);
    }
  else if(met==MOST_FREQUENT)
    {
      int i, j, numcath=0, maxnum=0;
      struct double_cathegory
      {
	double val;
	int num;
      };
      double_cathegory *cath=new double_cathegory[len];
      
      for(i=0;i<len;i++)
	cath[i].num=0;

      for(i=0;i<len;i++)
	{
	  int found=0;

	  for(j=0;j<numcath && !found;j++)
	    if(almost_equal(cath[j].val, data[i]))
	      {
		found=1;
		cath[j].num++;
	      }

	  if(!found)
	    {
	      cath[numcath].val=data[i];
	      cath[numcath].num++;
	      
	      numcath++;
	    }
	}
      
      for(j=0;j<numcath;j++)
	if(maxnum < cath[j].num)
	  {
	    val = cath[j].val;
	    maxnum = cath[j].num;
	  }

      delete [] cath;
    }

  return val;
}


// Returns the one-step auto-correlation of a sample 'x'
double get_auto_correlation(double *x, int len)
{
  double cov=0;
  double var=find_statistics(x,len,VARIATION);
  double mu=find_statistics(x,len,MEAN);

  for(int i=1;i<len;i++)
    cov+=(x[i-1]-mu)*(x[i]-mu);
  cov/=double(len-1);
  
  return cov/var;
}




// *  Description:       To perform a linear interpolation 
// *-----------------------------------------------------------------------*
// *  Programmer: Lars A. Roald, HD
// *  Revised by: Maitrayi Sabaratnam                      HydDate: 30.6.1993
// *  Revised by: Trond Reitan                             HydDate:  9.9.2016
// *  Changes:    Fortran  to C
// *-----------------------------------------------------------------------*
// *  In-parameters:                                                       *
// *  Ngame       Type     Descri1ption                                      *
// *  X0         real     Lowest X-value in the interval defining the      *
// *                      slope of the linear relation                     *
// *  Y0         real     Lowest Y-value                                   *
// *  X1         real     Highest X-value                                  *
// *  Y1         real     Highest Y-value                                  *
// *  X          real     X-value to be interpolated
// *************************************************************************
double linear(double x0, double y0, double x1, double y1, double x) 
{
  if (x1 - x0 == 0.0 || y0==MISSING_VALUE || y1==MISSING_VALUE) 
    return MISSING_VALUE;
  else 
    return (((y1-y0) / (x1-x0)) * (x-x0) + y0);
}



// *  Description:       To perform a linear interpolation 
// *-----------------------------------------------------------------------*
// *  Programmer: Lars A. Roald, HD
// *  Revised by: Maitrayi Sabaratnam                      HydDate: 30.6.1993
// *  Revised by: Trond Reitan                             HydDate:  9.9.2016
// *  Changes:    Fortran  to C
// *-----------------------------------------------------------------------*
// *  In-parameters:                                                       *
// *  Name       Type     Description                                      *
// *  X0         real     Lowest X-value in the interval defining the      *
// *                      slope of the linear relation                     *
// *  Y0         real     Lowest Y-value                                   *
// *  X1         real     Highest X-value                                  *
// *  Y1         real     Highest Y-value                                  *
// *  X          real     X-value to be interpolated
// *************************************************************************
double linear(HydDateTime x0, double y0, HydDateTime x1, double y1, HydDateTime x) 
{
  if ((x1 - x0) == 0 || y0 == MISSING_VALUE || y1 == MISSING_VALUE) 
    return MISSING_VALUE;
  else 
    return (((y1-y0) / double(x1-x0)) * double(x-x0) + y0);
}



int almost_equal(double v1, double v2)
{
  if(fabs(v1-v2)<1e-20)
    return 1;

  if(fabs(v1)>0.999999*fabs(v2) && fabs(v1)<1.0000001*fabs(v2) &&
     ((v1<0.0 && v2< 0.0) || (v1>0.0 && v2>0.0)))
    return 1;

  return 0;
}


//////////////////////////////////////////
// Matrix operations:
//////////////////////////////////////////

double **Make_matrix(int rows, int columns)
{
  int i,j;
  double **ret=new double*[rows];
  
  for(i=0;i<rows;i++)
    {
      ret[i]=new double[columns];
      for(j=0;j<columns;j++)
	ret[i][j]=0.0;
    }
  
  return ret;
}

Complex **Make_Complex_matrix(int rows, int columns)
{
  int i,j;
  Complex **ret=new Complex*[rows];
  
  for(i=0;i<rows;i++)
    {
      ret[i]=new Complex[columns];
      for(j=0;j<columns;j++)
	ret[i][j]=0.0;
    }

  return ret;
}

bool check_matrix(double **X, int n, int m)
{
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      if(!(X[i][j]>=-1e+200 && X[i][j]<=1e+200))
	return false;
  return true;
}

void show_mat(const char *name, double **X, int n, int m)
{
  int i,j;

#ifdef MAIN
  printf("%10s=(",name);
  for(i=0;i<n;i++)
    {
      for(j=0;j<m;j++)
	printf("%9.5f  ", X[i][j]);
      if(i<(n-1))
	printf(" )\n%11s("," ");
      else
	printf(" )\n\n");
    }
#else
  Rcout << name << "=(";
  for(i=0;i<n;i++)
    {
      for(j=0;j<m;j++)
	Rcout << X[i][j] << " ";
      if(i<(n-1))
	Rcout << ")" << std::endl << "           (";
      else
	Rcout << " )" << std::endl << std::endl;
    }
#endif // MAIN
}

void show_Complex_mat(const char *name, Complex **X, int n, int m)
{
  int i,j;

#ifdef MAIN
  printf("%10s=(",name);
  for(i=0;i<n;i++)
    {
      for(j=0;j<m;j++)
	printf("%9.5f+%9.5fi  ", X[i][j].Re(),X[i][j].Im());
      if(i<(n-1))
	printf(" )\n%11s("," ");
      else
	printf(" )\n\n");
    }
#else
  Rcout << name << "=(";
  for(i=0;i<n;i++)
    {
      for(j=0;j<m;j++)
	Rcout << X[i][j].Re() << "+" << X[i][j].Im() << "i ";
      if(i<(n-1))
	Rcout << ")" << std::endl << "           (";
      else
	Rcout << " )" << std::endl << std::endl;
    }
#endif // MAIN
}

void show_mat_R(const char *name, double **X, int n, int m)
{
  int i,j;
  
  char str[1000];
  snprintf(str,999,"%10s=(",name);
#ifdef MAIN
  cout << str;
#else
  Rcout << str;
#endif // MAIN
  for(i=0;i<n;i++)
    {
      for(j=0;j<m;j++)
	{
	  snprintf(str,999,"%9.5f  ", X[i][j]);
#ifdef MAIN
	  cout << str;
#else
	  Rcout << str;
#endif // MAIN
	}
      if(i<(n-1))
        snprintf(str,999," )\n%11s("," ");
      else
	snprintf(str,99," )\n\n");
#ifdef MAIN
      cout << str;
#else
      Rcout << str;
#endif // MAIN
    }
}

void show_Complex_mat_R(const char *name, Complex **X, int n, int m)
{
  int i,j;
  
  char str[1000];
  snprintf(str,999,"%10s=(",name);
#ifdef MAIN
  cout << str;
#else
  Rcout << str;
#endif // MAIN
  for(i=0;i<n;i++)
    {
      for(j=0;j<m;j++)
	{
	  snprintf(str,999,"%9.5f+%9.5fi  ", X[i][j].Re(),X[i][j].Im());
#ifdef MAIN
	  cout << str;
#else
	  Rcout << str;
#endif // MAIN
	}
      if(i<(n-1))
	snprintf(str,999," )\n%11s("," ");
      else
	snprintf(str,999," )\n\n");
#ifdef MAIN
      cout << str;
#else
      Rcout << str;
#endif // MAIN
    }
}

void show_vec(const char *name, double *X, int n)
{
  int j;
  
#ifdef MAIN
  printf("%10s=(",name);
  for(j=0;j<n;j++)
    printf("%9.5f  ", X[j]);
  printf(")\n");
#else
  Rcout << name << "=(";
  for(j=0;j<n;j++)
    Rcout << X[j] << "  ";
  Rcout << ")" << std::endl;
#endif // MAIN
}

void show_Complex_vec(const char *name, Complex *X, int n)
{
  int j;
  
#ifdef MAIN
  printf("%10s=(",name);
  for(j=0;j<n;j++)
    printf("%9.5f+%9.5fi  ", X[j].Re(),X[j].Im());
  printf(")\n");
#else
  Rcout << name << "=(";
  for(j=0;j<n;j++)
    Rcout << X[j].Re() << "+" << X[j].Im() << "i  ";
  Rcout << ")" << std::endl;
#endif // MAIN
}

#ifdef __cplusplus
extern "C" {
#endif

double **matrix_inverse(double **X, int dimension)
{
  //Rcout << "matrix_inverse" << std::endl;

#ifdef DETAILED_TIMERS
  timers[10][0]=clock();
#endif // DETAILED_TIMERS
  
  int i,j,n=dimension;
  double *lapack_matrix=new double[n*n];
  double *ret_matrix=new double[n*n];
  int info=0;
  double **ret=Make_matrix(n,n);
  int *ipv=new int[n];
  
  for(i=0;i<n;i++)
    ipv[i]=i+1;
  
  for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
	{
	  lapack_matrix[i+j*n]=X[i][j];
	  if(i==j)
	    ret_matrix[i+j*n]=1.0;
	  else
	    ret_matrix[i+j*n]=0.0;
	}
    }
  
  int n1=n,n2=n,n3=n,n4=n;

  /* Reminder:
  subroutine dgesv 	( 	integer  	N,
		integer  	NRHS,
		double precision, dimension( lda, * )  	A,
		integer  	LDA,
		integer, dimension( * )  	IPIV,
		double precision, dimension( ldb, * )  	B,
		integer  	LDB,
		integer  	INFO 
	)
	N	

  Where:
          N is INTEGER
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0.

[in]	NRHS	

          NRHS is INTEGER
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0.

[in,out]	A	

          A is DOUBLE PRECISION array, dimension (LDA,N)
          On entry, the N-by-N coefficient matrix A.
          On exit, the factors L and U from the factorization
          A = P*L*U; the unit diagonal elements of L are not stored.

[in]	LDA	

          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).

[out]	IPIV	

          IPIV is INTEGER array, dimension (N)
          The pivot indices that define the permutation matrix P;
          row i of the matrix was interchanged with row IPIV(i).

[in,out]	B	

          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
          On entry, the N-by-NRHS matrix of right hand side matrix B.
          On exit, if INFO = 0, the N-by-NRHS solution matrix X.

[in]	LDB	

          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).

[out]	INFO	

          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value
          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
                has been completed, but the factor U is exactly
                singular, so the solution could not be computed.
  */

  dgesv_(&n1,&n2,lapack_matrix,&n3,ipv,ret_matrix,&n4,&info);
  //LAPACKE_dgesv(LAPACK_COL_MAJOR, n1, n2, lapack_matrix, n3, ipv, info, n);
  
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      ret[i][j]=ret_matrix[i+j*n];
  
  delete [] lapack_matrix;
  delete [] ret_matrix;
  delete [] ipv;
  
#ifdef DETAILED_TIMERS
  timers[10][1]=clock();
  timers[10][2]+=(timers[10][1]-timers[10][0]);
#endif // DETAILED_TIMERS
  
  return ret;
}

#ifdef __cplusplus
}
#endif

Complex **Complex_matrix_inverse(Complex **X, int dimension)
{
  //Rcout << "Complex_matrix_inverse" << std::endl;

#ifdef DETAILED_TIMERS
  timers[11][0]=clock();
#endif // DETAILED_TIMERS
  
  int i,j,k,n=dimension;
  Complex c(1,45.0,false);
  Complex **cX=Make_Complex_matrix(n,n);
  Complex **invcX=Make_Complex_matrix(n,n);
  Complex **ret=Make_Complex_matrix(n,n);
  
  // Use the fact that inv(X)=c*inv(cX) to deal with degenerated real- and 
  // imaginary parts of the X matrix:

  //show_complex_mat("X",X,n,n);

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      cX[i][j]=c*X[i][j];

  //show_complex_mat("cX",cX,n,n);

  // Use that if X=A+Bi, where A and B are real matrixes, then
  // inv(X)=inv(A+B*inv(A)*B)-inv(B+A*inv(B)*A)*i
  
  double **A=Make_matrix(n,n),**B=Make_matrix(n,n);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      {
	A[i][j]=cX[i][j].Re();
	B[i][j]=cX[i][j].Im();
      }
  
  double **invB=matrix_inverse(B,n);
  //show_mat("B",B,n,n);
  //show_mat("invB",invB,n,n);
  double **invA=matrix_inverse(A,n);
  //show_mat("A",A,n,n);
  //show_mat("invA",invA,n,n);
  
  double **BinvA=Make_matrix(n,n), **BinvAB=Make_matrix(n,n);
  double **AinvB=Make_matrix(n,n), **AinvBA=Make_matrix(n,n);
  double **invRe=Make_matrix(n,n), **invIm=Make_matrix(n,n);

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      for(k=0;k<n;k++)
	{
	  BinvA[i][j]+=B[i][k]*invA[k][j];
	  AinvB[i][j]+=A[i][k]*invB[k][j];
	}

  //show_mat("BinvA",BinvA,n,n);
  //show_mat("AinvB",AinvB,n,n);

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)	
      for(k=0;k<n;k++)
	{
	  BinvAB[i][j]+=BinvA[i][k]*B[k][j];
	  AinvBA[i][j]+=AinvB[i][k]*A[k][j];
	}
  
  //show_mat("BinvAB",BinvA,n,n);
  //show_mat("AinvBA",AinvB,n,n);

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      {
	invRe[i][j]=A[i][j]+BinvAB[i][j];
	invIm[i][j]=B[i][j]+AinvBA[i][j];
      }
  
  //show_mat("invRe",invRe,n,n);
  //show_mat("invIm",invIm,n,n);


  double **plusRe=matrix_inverse(invRe,n), **minusIm=matrix_inverse(invIm,n);
  
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      {
	Complex newval(plusRe[i][j],-minusIm[i][j]);
	invcX[i][j]=newval;
      }
  
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      ret[i][j]=c*invcX[i][j];
  
  doubledelete(cX,n);
  doubledelete(invcX,n);
  doubledelete(A,n);
  doubledelete(B,n);
  doubledelete(invA,n);
  doubledelete(invB,n);
  doubledelete(BinvA,n);
  doubledelete(AinvB,n);
  doubledelete(BinvAB,n);
  doubledelete(AinvBA,n);
  doubledelete(invRe,n);
  doubledelete(invIm,n);
  doubledelete(plusRe,n);
  doubledelete(minusIm,n);

  
#ifdef DETAILED_TIMERS
  timers[11][1]=clock();
  timers[11][2]+=(timers[11][1]-timers[11][0]);
#endif // DETAILED_TIMERS
  
  
  return ret;
}


#ifdef __cplusplus
extern "C" {
#endif

  void call_dgeev(int n, double *X, double *wr, double *wi, double *vl, 
		  double *vr, double *work, int *info)
  {  
    int N=n,N2=n,N3=n,N4=n,N5=4*n;

    
    /*
      subroutine dgeev 	( 	character  	JOBVL,
		character  	JOBVR,
		integer  	N,
		double precision, dimension( lda, * )  	A,
		integer  	LDA,
		double precision, dimension( * )  	WR,
		double precision, dimension( * )  	WI,
		double precision, dimension( ldvl, * )  	VL,
		integer  	LDVL,
		double precision, dimension( ldvr, * )  	VR,
		integer  	LDVR,
		double precision, dimension( * )  	WORK,
		integer  	LWORK,
		integer  	INFO 
	)

[in]	JOBVL	

          JOBVL is CHARACTER*1
          = 'N': left eigenvectors of A are not computed;
          = 'V': left eigenvectors of A are computed.

[in]	JOBVR	

          JOBVR is CHARACTER*1
          = 'N': right eigenvectors of A are not computed;
          = 'V': right eigenvectors of A are computed.

[in]	N	

          N is INTEGER
          The order of the matrix A. N >= 0.

[in,out]	A	

          A is DOUBLE PRECISION array, dimension (LDA,N)
          On entry, the N-by-N matrix A.
          On exit, A has been overwritten.

[in]	LDA	

          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).

[out]	WR	

          WR is DOUBLE PRECISION array, dimension (N)

[out]	WI	

          WI is DOUBLE PRECISION array, dimension (N)
          WR and WI contain the real and imaginary parts,
          respectively, of the computed eigenvalues.  Complex
          conjugate pairs of eigenvalues appear consecutively
          with the eigenvalue having the positive imaginary part
          first.

[out]	VL	

          VL is DOUBLE PRECISION array, dimension (LDVL,N)
          If JOBVL = 'V', the left eigenvectors u(j) are stored one
          after another in the columns of VL, in the same order
          as their eigenvalues.
          If JOBVL = 'N', VL is not referenced.
          If the j-th eigenvalue is real, then u(j) = VL(:,j),
          the j-th column of VL.
          If the j-th and (j+1)-st eigenvalues form a complex
          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
          u(j+1) = VL(:,j) - i*VL(:,j+1).

[in]	LDVL	

          LDVL is INTEGER
          The leading dimension of the array VL.  LDVL >= 1; if
          JOBVL = 'V', LDVL >= N.

[out]	VR	

          VR is DOUBLE PRECISION array, dimension (LDVR,N)
          If JOBVR = 'V', the right eigenvectors v(j) are stored one
          after another in the columns of VR, in the same order
          as their eigenvalues.
          If JOBVR = 'N', VR is not referenced.
          If the j-th eigenvalue is real, then v(j) = VR(:,j),
          the j-th column of VR.
          If the j-th and (j+1)-st eigenvalues form a complex
          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
          v(j+1) = VR(:,j) - i*VR(:,j+1).

[in]	LDVR	

          LDVR is INTEGER
          The leading dimension of the array VR.  LDVR >= 1; if
          JOBVR = 'V', LDVR >= N.

[out]	WORK	

          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

[in]	LWORK	

          LWORK is INTEGER
          The dimension of the array WORK.  LWORK >= max(1,3*N), and
          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
          performance, LWORK must generally be larger.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.

[out]	INFO	

          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value.
          > 0:  if INFO = i, the QR algorithm failed to compute all the
                eigenvalues, and no eigenvectors have been computed;
                elements i+1:N of WR and WI contain eigenvalues which
                have converged.

    */
    
    //int irrel1=4*n,irrel2=4*n;
    dgeev_("V", "V", &N, X, &N2, wr, wi, vl, &N3, vr, &N4, work,
	   &N5, info, 1, 1);
			
    //LAPACKE_dgeev(LAPACK_COL_MAJOR, 'V', 'V', N, X, N2, wr, wi, vl, N3, vr, N4);  
  }

#ifdef __cplusplus
}
#endif


bool ml_started=false;
Complex **matrix_eigenvectors(double **A, int size, Complex **eigen_values,
			      Complex ***V_inv)
{
  int i,j,n=size;
  
#ifdef DETAILED_TIMERS
  timers[12][0]=clock();
#endif // DETAILED_TIMERS

  //Rcout << "matrix_eigenvectors" << std::endl;
  
  int info;
  double *vl=new double[n*n];
  double *vr=new double[n*n];
  double *wr=new double[n*n];
  double *wi=new double[n*n];
  double *work=new double[10*n*n];
  double *X=new double[n*n];
  
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      X[i+j*n]=A[i][j];
  
  call_dgeev(n, X, wr, wi, vl, vr, work, &info);

  double **vr_mat=Make_matrix(n,n);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      vr_mat[i][j]=vr[i+j*n];
  
  double **vr_re=Make_matrix(n,n), **vr_im=Make_matrix(n,n);
  
  for(i=0;i<n;i++) // traverse columns
    {
      double wi_abs=wi[i]>0.0 ? wi[i] : -wi[i];
      if(wi_abs>1e-9 && i<(n-1) && wi[i]==-wi[i+1])
	{
	  for(j=0;j<n;j++)
	    {
	      vr_re[j][i]=vr_re[j][i+1]=vr_mat[j][i];
	      vr_im[j][i]=vr_mat[j][i+1];
	      vr_im[j][i+1]=-vr_mat[j][i+1];
	    }
	  i++;
	}
      else
	{
	  for(j=0;j<n;j++)
	    {
	      vr_re[j][i]=vr_mat[j][i];
	      vr_im[j][i]=0.0;
	    }
	}
    }  

  Complex **vr_complex=Make_Complex_matrix(n,n);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      {
	Complex newvar(vr_re[i][j],vr_im[i][j]);
	vr_complex[i][j]=newvar;
      }
 
  Complex **vl_complex=Complex_matrix_inverse(vr_complex,n);
  
  Complex *lambda=new Complex[n];
  for(i=0;i<n;i++)
    {
      Complex newvar(wr[i],wi[i]);
      lambda[i]=newvar;
    }
  
  doubledelete(vr_mat,n);
  doubledelete(vr_re,n);
  doubledelete(vr_im,n);
  delete [] vl;
  delete [] vr;
  delete [] wr;
  delete [] wi;
  delete [] work;
  delete [] X;

  if(V_inv)
    *V_inv=vl_complex;
  else
    doubledelete(vl_complex,n);
  *eigen_values=lambda;
  
#ifdef DETAILED_TIMERS
  timers[12][1]=clock();
  timers[12][2]+=(timers[12][1]-timers[12][0]);
#endif // DETAILED_TIMERS

  return vr_complex;
  
}    

// check if a matrix has non-zero diagonal elements
bool Is_diagonal(Complex **A, int size)
{
  bool ret=true;
  
  for(int i=0;i<size && ret;i++)
    if(A[i][i].abs2()<1e-20)
      ret=false;
  
  return ret;
}

Complex *Complex_eigenvalues(double **A, int size)
{
  /*
  complex *val, **V=get_complex_eigenvector_matrix(A,size,&val,0,1000000);
  doubledelete(V,size);
  return val;
  */

  //Rcout << "Complex eigenvalues" << std::endl;
  
  Complex **V,**Vinv, *val;
  V=matrix_eigenvectors(A, size, &val, &Vinv);
  
  doubledelete(Vinv,size);
  doubledelete(V,size);
  return val;
}

double *double_eigenvalues(double **A, int size)
{
  //Rcout << "double_eigenvalues" << std::endl;

  double *ret=new double[size];
  Complex *val=Complex_eigenvalues(A,size);
  for(int i=0;i<size;i++)
    ret[i]=val[i].Re();
  
  delete [] val;
  return ret;
}

// find the determinant of a given matrix using eigenvector analysis;
double matrix_determinant(double **matrix_, int dimension)
{
  //Rcout << "matrix determinant" << std::endl;
  Complex *val=Complex_eigenvalues(matrix_,dimension);
  Complex det=1.0;
  
  for(int i=0;i<dimension;i++)
    det*=val[i];
  
  delete [] val;
  
  return det.Re();
}

Complex log_matrix_determinant(double **matrix_, int dimension)
{
  //Rcout << "log matrix determinant" << std::endl;
  Complex *val=Complex_eigenvalues(matrix_,dimension);
  Complex logdet=0.0;
  
  for(int i=0;i<dimension;i++)
    logdet+=complex_log(val[i]);
  
  delete [] val;
  
  return logdet;
}



double **get_cholesky(double **matrix_, int dimension)
{
  //Rcout << "get_cholesky" << std::endl;
  
  double **ret=Make_matrix(dimension,dimension);
  int i,j,n=dimension,lda=n,info;
  size_t size_of_something_mysterious=0;
  double *X=new double[n*n];
  
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      X[i+j*n]=matrix_[i][j];
  
  dpotrf_("U",&n,X,&lda,&info,size_of_something_mysterious);
  //LAPACKE_dpotrf(LAPACK_COL_MAJOR,'U',lapack_int(n),X,lda);
  
  for(i=0;i<n;i++)
    for(j=0;j<i;j++)
      X[i+j*n]=0.0;
  
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      ret[i][j]=X[i+j*n];
  
  delete [] X;
  return ret;
}



double get_random_gauss(void)
{
  return sqrt(-2.0*log(drand()))*sin(2.0*M_PI*drand());
}

double **multinormal_sample(int number_of_samples, double *expectation, 
			    double **sigma, int dimension, 
			    bool dont_keep_cholesky)
{
  // DEBUG static gsl_rng *rptr=gsl_rng_alloc(gsl_rng_rand48); 

  //Rcout << "multinormal_sample" << std::endl;
  
#ifdef DETAILED_TIMERS
  timers[13][0]=clock();
#endif // DETAILED_TIMERS

  if(dimension==1)
    {
      double **ret=new double*[number_of_samples];
      double s=sqrt(sigma[0][0]);
      
      for(int i=0;i<number_of_samples;i++)
	{
	  ret[i]=new double[1];
	  ret[i][0]=expectation[0]+s*get_random_gauss(); 
	  // DEBUG ret[i][0]=expectation[0]+s*gsl_ran_ugaussian(rptr);
	}
      
      return ret;
    }
  
  int i,j,k,len=dimension;
  
  if(dont_keep_cholesky)
    {
      double **c=get_cholesky(sigma, dimension);
      if(!c)
	return NULL;

      /* DEBUG
      double **c2=cholesky(sigma, dimension);
      
      for(i=0;i<len;i++)
	for(j=0;j<len;j++)
	  if(ABSVAL((c[i][j]-c2[i][j]))>0.001)
	    {
	      show_mat("sigma",sigma,len,len);
	      show_mat("c",c,len,len);
	      show_mat("c2",c2,len,len);
	    }
      doubledelete(c2,len);
      */

      double *u=new double[len], *x=new double[len];
      double **ret=new double*[number_of_samples];
      // fetch the samples;
      for(i=0;i<number_of_samples;i++)
	{
	  ret[i]=new double[len];
	  
	  // make a vector U~N(0,ident(len))
	  for(j=0;j<len;j++)
	    u[j]=get_random_gauss();  
	  // DEBUG u[j]=gsl_ran_gaussian(rptr, 1.0);
      
	  for(j=0;j<len;j++)
	    {
	      x[j]=0.0;
	      for(k=0;k<len;k++)
		x[j]+=c[k][j]*u[k];
	      // makes a vector with the correct correlation structure
	    }
	  
	  for(j=0;j<len;j++)
	    ret[i][j]=x[j]+expectation[j];
	  // put it into the return array
	}

      delete [] x;
      delete [] u;
      doubledelete(c,dimension);

      return ret;
    }

  static int visited=0;
  static double **sig=NULL;
  static int dim=0;
  static double **c=NULL;
  
  if(visited || sigma!=sig ||  dim!=dimension)
    {
      if(c)
	doubledelete(c, dim);
      c=get_cholesky(sigma, dimension); 
      if(!c)
	return NULL;

      /* DEBUG
      double **c2=cholesky(sigma, dimension);
      
      for(i=0;i<len;i++)
	for(j=0;j<len;j++)
	  if(ABSVAL((c[i][j]-c2[i][j]))>0.001)
	    {
	      show_mat("sigma",sigma,len,len);
	      show_mat("c",c,len,len);
	      show_mat("c2",c2,len,len);
	    }
      doubledelete(c2,len);
      */
      
      sig=sigma;
      dim=dimension;
    }

  double *u=new double[len], *x=new double[len];
  double **ret=new double*[number_of_samples];
  // fetch the samples;
  for(i=0;i<number_of_samples;i++)
    {
      ret[i]=new double[len];

      // make a vector U~N(0,ident(len))
      for(j=0;j<len;j++)
	u[j]=get_random_gauss(); 
      // DEBUG u[j]=gsl_ran_gaussian(rptr, 1.0);
      
      for(j=0;j<len;j++)
	{
	  x[j]=0.0;
	  for(k=0;k<len;k++)
	    x[j]+=c[k][j]*u[k];
	  // makes a vector with the correct correlation structure
	}

      for(j=0;j<len;j++)
	ret[i][j]=x[j]+expectation[j];
      // put it into the return array
    }

  delete [] x;
  delete [] u;
  visited=1;

#ifdef DETAILED_TIMERS
  timers[13][1]=clock();
  timers[13][2]+=(timers[13][1]-timers[13][0]);
#endif // DETAILED_TIMERS

  return ret;
}


double standard_normal_cdf(double x)
{
  return 0.5*(1.0+erf(x/sqrt(2)));
}

double normal_cdf(double x, double mean, double sigma)
{
  return standard_normal_cdf((x-mean)/sigma);
}


double chisq_deg1_cdf(double x)
{
  if(x<0.0)
    return MISSING_VALUE;
  
  return 2*standard_normal_cdf(sqrt(x))-1.0;
}

// Newtons method for a function with non-specified derivated.
// Calculates the derivated numerically
// Parameters; func  - the function
//             start - initializes the argument
//             enddiff - The iterations are stopped when 
//                       the diffierence from one iteration and the next
//                       is less than this value
//             maxiter - The maximal number of iterations. Returns the number 
//                       of iterations actually used
//             func_value - The value the function should have
//             steplength - The speed of the operation. In standard
//                          newton's method this value should be 1.0
//                          0.5 chosen for stability
double newtons_method(double (*func)(double), double start, double enddiff, 
		      int &maxiter, double func_value, double steplength, 
		      double minvalue,double maxvalue)
{
  double x=start, stepsize=steplength;
  int i=0;

  // Loop for as long as the number of iterations is less than the maximum,
  // and the function value is outside the interval we're interested in;
  while(i<maxiter && ABSVAL(((*func)(x) - func_value))>enddiff)
    {
      // Calculate the numeric derivative;
      double df = (((*func)(x+enddiff/20.0))-((*func)(x-enddiff/20.0)))/
	(enddiff/10.0);

      double step=(((*func)(x))-func_value)/df;
      
      if(minvalue!=MISSING_VALUE || maxvalue!=MISSING_VALUE)
	while((minvalue!=MISSING_VALUE && x-stepsize*step<=minvalue) ||
	      (maxvalue!=MISSING_VALUE && x-stepsize*step>=maxvalue))
	  stepsize/=2.0;
      // calculate the next argument from that;
      x = x - stepsize*step;
      stepsize=steplength;
      
      i++; // update the index
    }
  
  maxiter=i; // set the number of iterations used

  return x; // return the argument
}

double normal_invcdf(double p, double mean, double sigma)
{
  int maxiter=1000;
  double x=newtons_method(standard_normal_cdf, 0.0, 0.000001,
			  maxiter, p, 1.0, MISSING_VALUE, MISSING_VALUE);
  
  if(maxiter>=1000)
    return MISSING_VALUE;

  return mean+sigma*x;
}


double pdf_multinormal(double *x, double *expectation, double **sigma, int dimension,
		       bool logarithmic)
{
  //Rcout << "pdf_multinormal" << std::endl;

  
#ifdef DETAILED_TIMERS
  timers[14][0]=clock();
#endif // DETAILED_TIMERS
  
  
  double n=double(dimension);
  double ret=logarithmic ? -n/2.0*log(2.0*M_PI) : pow(2.0*M_PI, -n/2.0);
  double **inv_sigma=matrix_inverse(sigma, dimension);
  if(!inv_sigma)
    return MISSING_VALUE;
  
  int j,k,len=dimension;
  double SS=0.0;

  double det=logarithmic ? log_matrix_determinant(sigma,dimension).Re() : 
    matrix_determinant(sigma,dimension);

  if(det==MISSING_VALUE)
    {
      doubledelete(inv_sigma, dimension);
      return MISSING_VALUE;
    }
  
  for(j=0;j<len;j++)
    for(k=0;k<len;k++)
      SS += (x[j]-expectation[j])*inv_sigma[j][k]*(x[k]-expectation[k]);

  if(!logarithmic)
    ret*=exp(-0.5*SS);
  else
    ret-=0.5*SS;
  
  if(!logarithmic)
    ret/=sqrt(det);
  else
    ret-=0.5*det;

  doubledelete(inv_sigma, dimension);

#ifdef DETAILED_TIMERS
  timers[14][1]=clock();
  timers[14][2]+=(timers[14][1]-timers[14][0]);
#endif // DETAILED_TIMERS
  
  return ret;
}

double pdf_multinormal(double *x, double *expectation, 
		       int dimension, double **inv_sigma, double det_sigma,
		       bool logarithmic)
{
#ifdef DETAILED_TIMERS
  timers[14][0]=clock();
#endif // DETAILED_TIMERS
  
  double n=double(dimension);
  double ret=logarithmic ? -n/2.0*log(2.0*M_PI) : pow(2.0*M_PI, -n/2.0);
  int j,k,len=dimension;
  double SS=0.0;

  for(j=0;j<len;j++)
    for(k=0;k<len;k++)
      SS += (x[j]-expectation[j])*inv_sigma[j][k]*(x[k]-expectation[k]);

  if(!logarithmic)
    ret*=exp(-0.5*SS);
  else
    ret-=0.5*SS;
  
  if(!logarithmic)
    ret/=sqrt(det_sigma);
  else
    ret-=0.5*det_sigma; // it is assumed that the determinant is now 
			// given in logarithmic form

#ifdef DETAILED_TIMERS
  timers[14][1]=clock();
  timers[14][2]+=(timers[14][1]-timers[14][0]);
#endif // DETAILED_TIMERS
  
  return ret;
}


double pdf_multi_t(double *x, double *expectation, double **sigma, double nu,
		   int dimension,
		   bool logarithmic)
{
#ifdef DETAILED_TIMERS
  timers[14][0]=clock();
#endif // DETAILED_TIMERS
  
  double n=double(dimension);
  double ret=lgamma(0.5*(nu+n))-lgamma(0.5*nu)-0.5*n*log(M_PI*nu);
  if(!logarithmic)
    ret=exp(ret);

  //Rcout << "pdf_multi_t" << std::endl;
  
  double **inv_sigma=matrix_inverse(sigma, dimension);
  int j,k,len=dimension;
  double det=logarithmic ? log_matrix_determinant(sigma,dimension).Re() : 
    matrix_determinant(sigma,dimension);
  double SS=0.0;

  for(j=0;j<len;j++)
    for(k=0;k<len;k++)
      SS += (x[j]-expectation[j])*inv_sigma[j][k]*(x[k]-expectation[k]);

  if(!logarithmic)
    ret*=pow(1.0+SS/nu,-0.5*(nu+n));
  else
    ret-=0.5*(nu+n)*log(1.0+SS/nu);
  
  if(!logarithmic)
    ret/=sqrt(det);
  else
    ret-=0.5*det;

  doubledelete(inv_sigma, dimension);

#ifdef DETAILED_TIMERS
  timers[14][1]=clock();
  timers[14][2]+=(timers[14][1]-timers[14][0]);
#endif // DETAILED_TIMERS
  
  return ret;
}

double pdf_multi_t(double *x, double *expectation, double nu,
		   int dimension, double **inv_sigma, double det_sigma,
		   bool logarithmic)
{
#ifdef DETAILED_TIMERS
  timers[14][0]=clock();
#endif // DETAILED_TIMERS
  
  double n=double(dimension);
  double ret=lgamma(0.5*(nu+n))-lgamma(0.5*nu)-0.5*n*log(M_PI*nu);
  if(!logarithmic)
    ret=exp(ret);
  int j,k,len=dimension;
  double SS=0.0;

  for(j=0;j<len;j++)
    for(k=0;k<len;k++)
      SS += (x[j]-expectation[j])*inv_sigma[j][k]*(x[k]-expectation[k]);

  if(!logarithmic)
    ret*=pow(1.0+SS/nu,-0.5*(nu+n));
  else
    ret-=0.5*(nu+n)*log(1.0+SS/nu);
  
  if(!logarithmic)
    ret/=sqrt(det_sigma);
  else
    ret-=0.5*det_sigma; // it is assumed that the determinant is now 
			// given in logarithmic form

#ifdef DETAILED_TIMERS
  timers[14][1]=clock();
  timers[14][2]+=(timers[14][1]-timers[14][0]);
#endif // DETAILED_TIMERS
  
  return ret;
}

double *mean_of_vectors(int number_of_samples, double **samples, int dimension)
{
#ifdef DETAILED_TIMERS
  timers[15][0]=clock();
#endif // DETAILED_TIMERS
  
  int i,j,len=dimension;
  double *ret=new double[len];

  for(j=0;j<len;j++)
    ret[j]=0.0;

  // Make the vector of sum(samples);
  for(i=0;i<number_of_samples;i++)
    for(j=0;j<len;j++)
      ret[j]+=samples[i][j];
  
  // Divide by the number of samples;
  for(j=0;j<len;j++)
    ret[j]/=double(number_of_samples);
  
#ifdef DETAILED_TIMERS
  timers[15][1]=clock();
  timers[15][2]+=(timers[15][1]-timers[15][0]);
#endif // DETAILED_TIMERS
  
  return ret;
}


double **estimated_variance(int number_of_samples,
				 double **samples /* a sample of vectors */,
				 double *mean, int dimension)
{
#ifdef DETAILED_TIMERS
  timers[15][0]=clock();
#endif // DETAILED_TIMERS
  
  int i,j,k,len=dimension;
  double **sigma=new double*[len];

  for(j=0;j<len;j++)
    {
      sigma[j]=new double[len];
      for(k=0;k<len;k++)
	sigma[j][k]=0.0;
    }

  // fill out the upper triangle of sigma with transpose(x-mu)*(x-mu);
  for(i=0;i<number_of_samples;i++)
    for(j=0;j<len;j++)
      for(k=j;k<len;k++)
	sigma[j][k]+=(samples[i][j]-mean[j])*
	  (samples[i][k]-mean[k]);

  // Divide by the number of samples
  for(j=0;j<len;j++)
    for(k=j;k<len;k++)
      sigma[j][k]/=double(number_of_samples);

  // fill out the lower triangle of sigma to make it symmetric
  for(j=0;j<len;j++)
    for(k=0;k<j;k++)
      sigma[j][k]=sigma[k][j];

#ifdef DETAILED_TIMERS
  timers[15][1]=clock();
  timers[15][2]+=(timers[15][1]-timers[15][0]);
#endif // DETAILED_TIMERS
  
  return sigma;
}


double **estimated_variance(int number_of_samples,
				 double **samples /* a sample of vectors */,
				 int dimension)
{
  double *mean=mean_of_vectors(number_of_samples, samples, dimension);
  double **ret=estimated_variance(number_of_samples, samples, mean, dimension);

  delete [] mean;

  return ret;
}





// ****************************************************
// Quasi-Newton optimizer taken from hydrasub-library
// ****************************************************


double *quasi_newton(double (*func)(double *), 
		     double *(*derivated)(double *), 
		     int dimension,
		     double *start, double enddiff, 
		     int &maxiter, bool do_min=true,
		     bool hillclimb=false,
		     bool silent=true);
double *quasi_newton(double (*func)(double *), 
		     int dimension,
		     double *start, double enddiff, 
		     int &maxiter, bool do_min=true,
		     bool hillclimb=false,
		     bool silent=true);

// Static global variable used inside the multivariate
// Newton-Raphson methods;
double (*quasi_func)(double *);
static int quasi_dimension;

// derivated of the last function used by multinewtonraphson
double *quasi_numeric_derivated(double *x)
{
  int i,j;
  double fval=(*quasi_func)(x); // fetch the function value
  double *df=new double[quasi_dimension]; // derivated vector
  double *x2=new double[quasi_dimension]; // argument vector buffer
  double delta=1e-7; // precision

  // traverse the dimensions
  for(i=0;i<quasi_dimension;i++)
    {
      // Set the argument vecotr, going a little to the left for one index, i;

      for(j=0;j<quasi_dimension;j++)
	if(i!=j)
	  x2[j]=x[j];
	else
	  x2[j]=x[j]+delta;

      // Calculate the derivated as the scaled difference between the
      // function value before and after the step;
      df[i]=((*quasi_func)(x2)-fval)/delta;
    }

  // cleanup;
  delete [] x2;

  return df; // return the derivated vector
}


double *quasi_newton(double (*func)(double *), 
		     double *(*derivated)(double *), 
		     int dimension,
		     double *start, double enddiff, 
		     int &maxiter, bool do_min,
		     bool hillclimb, bool silent)
{
#ifdef DETAILED_TIMERS
  timers[23][0]=clock();
#endif // DETAILED_TIMERS
  
  double *x=new double[dimension]; // argument vector
  double *x1=new double[dimension]; // argument vector buffer
  int i=0,j,k,l, m;

  /*
  if(!silent)
#ifdef MAIN
    cout << "quasi_newton, enddiff=" << enddiff << " maxiter=" <<
      maxiter << " do_min=" << do_min << " hillclimb=" << hillclimb <<
      std::endl;
#else
    Rcout << "quasi_newton, enddiff=" << enddiff << " maxiter=" <<
      maxiter << " do_min=" << do_min << " hillclimb=" << hillclimb <<
      std::endl;
#endif MAIN
  */
  
  for(k=0;k<dimension;k++)
    x[k]=start[k];

  // find the derivated vector;
  double f0=(*func)(x), f1=MISSING_VALUE;
  double *df=(*derivated)(x);
  double **H=new double*[dimension];
  double **Hnew=new double*[dimension];
  double *s=new double[dimension];
  double *p=new double[dimension];
  double *q=new double[dimension];
  double qHq,ptq,*Hq=new double[dimension],**pqt=new double*[dimension];

  // Set H=unity matrix
  for(k=0;k<dimension;k++)
    {
      H[k]=new double[dimension];
      pqt[k]=new double[dimension];
      Hnew[k]=new double[dimension];
      for(l=0;l<dimension;l++)
	if(k==l)
	  H[k][l]=1.0;
	else
	  H[k][l]=0.0;
    }

  double lambda[13]={1e-12, 1e-6, 0.001,
		     0.003,0.01,0.03,0.1,0.3,1.0,3.0,10.0,30.0,100.0};


  // Loop for as long as the number of iterations is less than the maximum,
  // and the function value is outside the interval we're interested in;
  while(i<maxiter && ABSVAL((f0-f1))>enddiff)
    { 
      df=(*derivated)(x); // fetch the derivated vector
      f0=(*func)(x);

      
      // s=H*df;
      if(!hillclimb)
	{
	  for(k=0;k<dimension;k++)
	    {
	      s[k]=0.0;
	      for(l=0;l<dimension;l++)
		s[k]+=H[k][l]*df[l];
	    }
	  if(!do_min)
	    for(k=0;k<dimension;k++)
	      s[k]=-s[k];
	}
      else // hillclimbing
	{
	  for(k=0;k<dimension;k++)
	    s[k]=df[k];
	  double norm=0;
	  for(k=0;k<dimension;k++)
	    norm+=s[k]*s[k];
	  norm=sqrt(norm);
	  for(k=0;k<dimension;k++)
	    s[k]/=norm;
	}

      //cout << H[0][0] << endl;
      //cout << s[0] << endl;

      // Minimize/maximize func(x-lambda*s)
      int mindex=0;
      double llmax=0.0;
      for(k=0;k<dimension;k++)
	x1[k]=x[k]-lambda[0]*s[k];
      double ex=(*func)(x1); 
      double sign=1.0;
      for(j=1;j<13;j++)
	{
	  for(k=0;k<dimension;k++)
	    x1[k]=x[k]-lambda[j]*s[k];
	  double f2=(*func)(x1);
	  if(f2!=MISSING_VALUE &&
	     (ex==MISSING_VALUE  || (do_min && f2<ex) || (!do_min && f2>ex)))
	    {
	      sign=1.0;
	      mindex=j;
	      ex=f2;
	      llmax=lambda[j];
	    }
	}
      if(mindex==12)
	{
	  double ll=lambda[12];
	  double f2=ex;
	  bool docont=true;
	  
	  do
	    {
	      docont=false;
	      ll*=10.0;
	      for(k=0;k<dimension;k++)
		x1[k]=x[k]-ll*s[k];
	      f2=(*func)(x1);
	      if(f2!=MISSING_VALUE &&
		 ((do_min && f2<ex) || (!do_min && f2>ex)))
		{
		  sign=1.0;
		  ex=f2;
		  llmax=ll;
		  docont=true;
		}
	    } while(docont);
	}
      
      for(j=0;j<13;j++)
	{
	  for(k=0;k<dimension;k++)
	    x1[k]=x[k]+lambda[j]*s[k];
	  double f2=(*func)(x1);
	  if(f2!=MISSING_VALUE && 
	     (ex==MISSING_VALUE  || (do_min && f2<ex) || (!do_min && f2>ex)))
	    {
	      sign=-1.0;
	      mindex=j;
	      ex=f2;
	      llmax=lambda[j];
	    }
	}
      if(sign<0.0 && mindex==12)
	{
	  double ll=lambda[12];
	  double f2=ex;
	  bool docont=true;
	  
	  do
	    {
	      docont=false;
	      ll*=10.0;
	      for(k=0;k<dimension;k++)
		x1[k]=x[k]+ll*s[k];
	      f2=(*func)(x1);
	      if(f2!=MISSING_VALUE &&
		 ((do_min && f2<ex) || (!do_min && f2>ex)))
		{
		  sign=-1.0;
		  ex=f2;
		  llmax=ll;
		  docont=true;
		}
	    } while(docont);
	}
      
      for(k=0;k<dimension;k++)
	x1[k]=x[k]-sign*llmax*s[k];
      double *df2=(*derivated)(x1); 

      //cout << lambda[mindex] << " " << x1[0] << endl;

      qHq=0.0; 
      ptq=0.0;
      for(k=0;k<dimension;k++)
	{
	  p[k]=x1[k]-x[k];
	  q[k]=df2[k]-df[k];
	}

      for(k=0;k<dimension;k++)
	{
	  Hq[k]=0.0;
	  for(l=0;l<dimension;l++)
	    {
	      Hq[k]+=H[k][l]*q[l];
	      pqt[k][l]=p[k]*q[l];
	    }
	  qHq+=q[k]*Hq[k];
	  ptq+=p[k]*q[k];
	}

      for(k=0;k<dimension;k++)
	for(l=0;l<dimension;l++)
	  Hnew[k][l]=H[k][l];

      double fac1=(1.0+qHq/ptq)/ptq;
      for(k=0;k<dimension;k++)
	for(l=0;l<dimension;l++)
	  Hnew[k][l]+=fac1*p[k]*p[l];
      
      for(k=0;k<dimension;k++)
	for(l=0;l<dimension;l++)
	  for(m=0;m<dimension;m++)
	    Hnew[k][l]-=(H[k][m]*pqt[m][l]+H[l][m]*pqt[m][k])/ptq;

      for(k=0;k<dimension;k++)
	{
	  x[k]=x1[k];
	  for(l=0;l<dimension;l++)
	    H[k][l]=Hnew[k][l];
	}
      
      f1=(*func)(x1);

      /*
      if(!silent)
	{
#ifdef MAIN
	  cout << "iteration " << i << " df0=" << df[0] << " f0=" << f0 <<
	    " f1=" << f1 << std::endl;
#else
	  Rcout << "iteration " << i << " df0=" << df[0] << " f0=" << f0 <<
	    " f1=" << f1 << std::endl;
#endif // MAIN
	  show_vec("x",x,dimension);
	  show_vec("x1",x1,dimension);
	  show_vec("s",s,dimension);
	}
      */
      
      delete [] df2;
      delete [] df;
      i++;
    }
  
  // cleanup;
  delete [] x1;
  doubledelete(H,dimension);
  doubledelete(Hnew,dimension);
  delete [] s;
  delete [] p;
  delete [] q;
  delete [] Hq;
  doubledelete(pqt,dimension);

  maxiter=i; // set the number of ieterations done
  
#ifdef DETAILED_TIMERS
  timers[23][1]=clock();
  timers[23][2]+=(timers[23][1]-timers[23][0]);
#endif // DETAILED_TIMERS
  
  return x; // return the argument vector  
}

double *quasi_newton(double (*func)(double *), 
		     int dimension,
		     double *start, double enddiff, 
		     int &maxiter, bool do_min,
		     bool hillclimb, bool silent)
{
  // Set the global function used for calculating numeric derivated to
  // the one given;
  quasi_func=func;
  quasi_dimension=dimension;

  // Use this derivated to calculate the Newton Raphson argument;
  return quasi_newton(func, quasi_numeric_derivated, dimension,
		      start, enddiff, maxiter, do_min,
		      hillclimb, silent);
}


// ****************************************************
// General-purpose MCMC, for use in residual analysis.
// Taken from hydrasub-library, file mcmc.C/H
// ****************************************************

double log_prob(double **data, unsigned int numrows, unsigned int numcolumns,
		unsigned int numparams, double *params,
		unsigned int num_hyperparameters, double *hyper_parameters,
		
		double (*log_prior)(unsigned int /*numparams*/, 
				    double* /* params */,
				    unsigned  int /* num_hyperparameters */,
				    double* /* hyper_parameters */),
		 
		double (*log_lik)(double** /*data*/, 
				  unsigned int /* numrows */, 
				  unsigned int /* numcolumns */,
				  unsigned int /*numparams*/, 
				  double* /* params */), 

		double T=1.0);

void sample_once(double **data, unsigned int numrows, unsigned int numcolumns,
		 unsigned int numparams, unsigned int numtemp, double **params, 
		 int *acc, double *rw, double *logprob,
		 double *T, int *swaps, double tempering_prob,
		 unsigned int num_hyperparameters, double *hyper_parameters,
		 
		 double (*log_prior)(unsigned int /*numparams*/, 
				     double* /* params */,
				     unsigned int /* num_hyperparameters */,
				     double* /* hyper_parameters */),
		 
		 double (*log_lik)(double** /*data*/, 
				   unsigned int /* numrows */, 
				   unsigned int /* numcolumns */,
				   unsigned int /*numparams*/, 
				   double* /* params */),
		 bool update_prev_log_prob=false ,
		 bool silent=true);
		


double **general_mcmc(unsigned int numsamples, unsigned int burnin,
		      unsigned int indep, unsigned int numtemp,

		      double **data,
		      unsigned int numrows, unsigned int numcolumns,
		      
		      unsigned int numparams, char **param_names, 
		      
		      double *T, double tempering_prob, 
		      
		      unsigned int num_hyperparameters,
		      double *hyper_parameters,
		      
		      void (*init)(unsigned int /* numparams */,
				   double* /* input parameter array*/, 
				   unsigned int /* num_hyperparameters */,
				   double* /* hyperparameters */), // init mechanism
		      // may possibly use the prior via hyper-parameters
		      
		      double (*log_prior)(unsigned int /*numparams*/, 
					  double* /* params */,
					  unsigned int /* num_hyperparameters */,
				 double* /* hyper_parameters */),
		      
		      double (*log_lik)(double** /*data*/, 
					unsigned int /* numrows */, 
					unsigned int /* numcolumns */,
					unsigned int /*numparams*/, 
					double* /* params */),
		      
		      bool silent=true,
		      bool update_prev_log_prob=false,
		      bool show_graphs=false);


// Importance sampling method for calculating model likelihood
// using a multinormal approximation of the posterior samples
// as the proposal distribution.
double log_model_likelihood_multinormal(unsigned int num_importance_samples,
			    
					double **mcmc_params, 
					unsigned int numparams,
					unsigned int numsamples,
					
					double **data,
					unsigned int numrows,
					unsigned int numcolumns,
					
					unsigned int num_hyperparameters,
					double *hyper_parameters,
					
					double (*log_prior)(unsigned int /*numparams*/, 
							    double* /* params */,
							    unsigned int /* num_hyperparameters */,
							    double* /* hyper_parameters */
							    ),
					
					double (*log_lik)(double** /*data*/, 
							  unsigned int /* numrows */, 
							  unsigned int /* numcolumns */,
							  unsigned int /*numparams*/, 
							  double* /* params */),
					
			                bool silent=true, 
					bool dont_keep_cholesky=false);


// show_parameter_value : prints current parameter array to screen:
void show_parameter_value(double *param, unsigned int numparam, char **parnames,
			  int max_params=(int) MISSING_VALUE);


// show_mcmc_parameters: Uses the external programs
// 'vvgraph' and 'histogramme' from the hydrasub
// package to show the sampling for a given
// parameter. If these are not present, this
// methods should be set on silent (only spacing
// between independent samples are shown).
// par: parameter samples
// N: number of samples
// parname: name of parameter
// silent: if set, turns on silent modus
void show_mcmc_parameter(double *par, unsigned int N, char *parname, 
			 bool silent=false, 
			 char *filestart=NULL);


#ifdef MAIN
// show_scatter: Shows scatterplots of the samples for
// two parameter.
// par1: parameter samples for parameter 1
// par2: parameter samples for parameter 2
// N: number of samples
// parname1: parameter name 1
// parname2: parameter name 2
void show_mcmc_scatter(double *par1, double *par2, int N, 
		       char *parname1, char *parname2, 
		       char *filestart=NULL);
#endif // MAIN

double log_prob(double **data, unsigned int numrows, unsigned int numcolumns,
		unsigned int numparams, double *params,
		unsigned int num_hyperparameters, double *hyper_parameters,
		
		double (*log_prior)(unsigned int /*numparams*/, 
				    double* /* params */,
				    unsigned int /* num_hyperparameters */,
				    double* /* hyper_parameters */),
		 
		double (*log_lik)(double** /*data*/, 
				  unsigned int /* numrows */, 
				  unsigned int /* numcolumns */,
				  unsigned int /*numparams*/, 
				  double* /* params */), 

		double T )
{
  double lp=log_prior(numparams,params,num_hyperparameters,hyper_parameters);

  if(!(lp>-1e+200 && lp<1e+200))
    return(-1e+200);

  double ll=log_lik(data,numrows,numcolumns,numparams,params);

  if(!(ll>-1e+200 && ll<1e+200))
    return(-1e+200);
  
  double pp=ll+lp;

  return pp/T;
}

void sample_once(double **data, unsigned int numrows, unsigned int numcolumns,
		 unsigned int numparams, unsigned int numtemp, double **params, 
		 int *acc, double *rw, double *logprob,
		 double *T, int *swaps, double tempering_prob,
		 unsigned int num_hyperparameters, double *hyper_parameters,
		 
		 double (*log_prior)(unsigned int /*numparams*/, 
				     double* /* params */,
				     unsigned int /* num_hyperparameters */,
				     double* /* hyper_parameters */),
		 
		 double (*log_lik)(double** /*data*/, 
				   unsigned int /* numrows */, 
				   unsigned int /* numcolumns */,
				   unsigned int /*numparams*/, 
				   double* /* params */),
		 bool update_prev_log_prob,
		 bool silent)
{
  unsigned int i;
  
  if(numtemp>1 && drand()<tempering_prob) // tempering swap?
    {
      // choose a random chain below the maximum:
      int index=(int) floor(double(numtemp-1)*drand());
      // Find the logprob when a switch takes place:
      double p1=log_prob(data,numrows,numcolumns,numparams,params[index],
			 num_hyperparameters, hyper_parameters,
			 log_prior, log_lik, T[index+1]);
      
      double p2=log_prob(data,numrows,numcolumns,numparams,params[index+1],
			 num_hyperparameters, hyper_parameters,
			 log_prior, log_lik, T[index]);
      
      if(update_prev_log_prob)
	{
	  logprob[index]=log_prob(data,numrows,numcolumns,numparams,
				  params[index],num_hyperparameters, 
				  hyper_parameters,
				  log_prior, log_lik, T[index]);
	  logprob[index+1]=log_prob(data,numrows,numcolumns,numparams,
				  params[index+1],num_hyperparameters, 
				  hyper_parameters,
				  log_prior, log_lik, T[index+1]);
	}

      // Acceptance-rejection term:
      if(log(drand())<p1+p2-logprob[index]-logprob[index+1])
	{
	  // Swap chains:
	  logprob[index+1]=p1;
	  logprob[index]=p2;
	  
	  // Switch the contents of the chains:
	  double *buffer=new double[numparams];
	  for(i=0;i<numparams;i++)
	    buffer[i]=params[index][i];
	  for(i=0;i<numparams;i++)
	    params[index][i]=params[index+1][i];
	  for(i=0;i<numparams;i++)
	    params[index+1][i]=buffer[i];
	  
	  swaps[index]++; // update the number of swaps
	  
	  if(!silent)
#ifdef MAIN
	    cout << " Swap chain " << index << " <-> chain " << index+1 << endl; 
#else
	    Rcout << " Swap chain " << index << " <-> chain " << index+1 << std::endl; 
#endif // MAIN	  

	  // cleanup
	  delete [] buffer;
	}
    }
  else // normal RW Metropolis sampling:
    for(unsigned int t=0;t<numtemp;t++) // traverse the tempering chains
      {
	double *newparams=new double[numparams];
	for(i=0;i<numparams;i++)
	  newparams[i]=params[t][i];
	
	double new_logprob;
	double sT=sqrt(T[t]); // tempering contribution to the 
	// random walk standard deviation
	
	for(i=0;i<numparams;i++)
	  {
	    newparams[i]+=rw[i]*sT*get_random_gauss();

	    // Update the logprob (log(prior*likelihood)):
	    new_logprob=log_prob(data,numrows,numcolumns,numparams,newparams,
				 num_hyperparameters, hyper_parameters,
				 log_prior, log_lik, T[t]);
	    
	    if(update_prev_log_prob)
	      logprob[t]=log_prob(data,numrows,numcolumns,numparams,params[t],
				  num_hyperparameters, hyper_parameters,
				  log_prior, log_lik, T[t]);

	    // Acceptance/rejection term:
	    if(log(drand())<new_logprob-logprob[t])
	      {
		// copy the new value over to the 'sample' structure:
		params[t][i]=newparams[i];
		// update the logprob:
		logprob[t]=new_logprob;
		// if this is the lowest chain (the chain of the
		// model we are interested in), then set the
		// acceptances indicator for this parameter:
		if(t==0)
		  acc[i]=1;
	      }
	    else // rejection
	      {
		// if this is the lowest chain (the chain of the
		// model we are interested in), then nullify the
		// acceptances indicator for this parameter:
		if(t==0)
		  acc[i]=0;
		// put the old sample value into the 
		// new sample structure
		newparams[i]=params[t][i];
	      }
	  }

	delete [] newparams;
      }
}


double **general_mcmc(unsigned int numsamples, unsigned int burnin,
		      unsigned int indep, unsigned int numtemp,

	      double **data, unsigned int numrows, unsigned int numcolumns,

	      unsigned int numparams, char **param_names, 

	      double *T, double tempering_prob, 

	      unsigned int num_hyperparameters, double *hyper_parameters,

	      void (*init)(unsigned int /* numparams */,
			   double* /* input parameter array*/, 
			   unsigned int /* num_hyperparameters */,
			   double* /* hyperparameters */), // init mechanism
	      // may possibly use the prior via hyper-parameters

	      double (*log_prior)(unsigned int /*numparams*/, 
				 double* /* params */,
				 unsigned int /* num_hyperparameters */,
				 double* /* hyper_parameters */),
	      
	      double (*log_lik)(double** /*data*/, 
				unsigned int /* numrows */, 
				unsigned int /* numcolumns */,
				unsigned int /*numparams*/, 
				double* /* params */),

	      bool silent,
	      bool update_prev_log_prob,
	      bool show_graphs)
{
  if(T[0]!=1.0)
    {
#ifdef MAIN
      cerr << "Usage in routine 'mcmc': first tempering temperature "
	"*must* be equal to 1!" << endl;
#else
      Rcout << "Usage in routine 'mcmc': first tempering temperature "
	"*must* be equal to 1!" << std::endl;
#endif // MAIN
      return NULL;
    }
  
  // Return matrix:
  double **ret=Make_matrix(numsamples,numparams);
  
  // indexes:
  unsigned int i,j,k,t;
  
  // Current set of samples (one for each tempering chain) and
  // acceptance indicators and random walk standard deviations:
  double **params=Make_matrix(numtemp,numparams);
  double *rw=new double[numparams];
  int *totalacc=new int[numparams];
  int *acc=new int[numparams];
  
  for(k=0;k<numparams;k++)
    {
      rw[k]=0.1;
      totalacc[k]=0;
    }
  
  // logprob, i.e. log(prior*likelihood):
  double *logprob=new double[numtemp];
  
  // Number of swaps between one chain and the next:
  int *swaps=new int[numtemp];
  for(t=0;t<numtemp;t++)
    swaps[t]=0;
  
#ifdef DETAILED_TIMERS
  timers[20][0]=clock();
#endif // DETAILED_TIMERS
  
  // **************************************************************
  // initialization of the samples:
  for(t=0;t<numtemp;t++) // traverse the chains
    {
      do
	{
	  init(numparams, params[t], num_hyperparameters, hyper_parameters);
	  
	  logprob[t]=log_prob(data, numrows, numcolumns,
			      numparams, params[t],
			      num_hyperparameters, hyper_parameters,
			      log_prior, log_lik, T[t]); 
	} while(!(logprob[t]>-1e+200 && logprob[t]<1e+200));
    }
  if(!silent)
#ifdef MAIN
    cout << "done getting initial values" << endl;
#else
    Rcout << "done getting initial values" << std::endl;
#endif // MAIN
    
#ifdef DETAILED_TIMERS
  timers[20][1]=clock();
  timers[20][2]+=(timers[20][1]-timers[20][0]);
#endif // DETAILED_TIMERS

  
#ifdef DETAILED_TIMERS
  timers[21][0]=clock();
#endif // DETAILED_TIMERS
  
  // **************************************************************
  // burn-in phase:
  if(burnin>0)
    {
      for(j=0;j<2;j++) // two sub-phases, with and without 
	// adaptation, done twice
	{
	  // sample without adaptation:
	  for(i=1;i<=burnin/4;i++)
	    {
	      sample_once(data,numrows,numcolumns,
			  numparams, numtemp, params,
			  acc, rw, logprob, T, swaps, tempering_prob,
			  num_hyperparameters, hyper_parameters,
			  log_prior, log_lik, update_prev_log_prob, silent);
	      if(!silent)
		{
#ifdef MAIN
		  cout << "burnin, j=" << j*2+1 << " of 4 i=" << i << 
		    " of " << burnin/4 << " logprob=" << logprob[0] << " ";
#else
		  Rcout << "burnin, j=" << j*2+1 << " of 4 i=" << i << 
		    " of " << burnin/4 << " logprob=" << logprob[0] << " ";
#endif // MAIN
		  show_parameter_value(params[0], numparams, param_names,10);
		}
	    }
	  
	  // update the acceptances:
	  for(k=0;k<numparams;k++)
	    totalacc[k]=0;
	  
	  // sample with adaptation:
	  for(i=1;i<=burnin/4;i++)
	    {
	      sample_once(data,numrows,numcolumns,
			  numparams, numtemp, params,
			  acc, rw, logprob, T, swaps, tempering_prob,
			  num_hyperparameters, hyper_parameters,
			  log_prior, log_lik, update_prev_log_prob, silent);
	      if(!silent)
		{
#ifdef MAIN
		  cout << "burnin, j=" << j*2+2 << " of 4 i=" << i << 
		    " of " << burnin/4 << " logprob=" << logprob[0] << " ";
#else
		  Rcout << "burnin, j=" << j*2+2 << " of 4 i=" << i << 
		    " of " << burnin/4 << " logprob=" << logprob[0] << " ";
#endif // MAIN
		  show_parameter_value(params[0], numparams, param_names,10);
		}

	      // update the total amount of acceptances for each parameter:
	      for(k=0;k<numparams;k++)
		totalacc[k]+=acc[k];
	      
	      // Adapt according to the acceptance rate each 100th sample:
	      if(i%100==0)
		{
		  // traverse the parameters:
		  for(k=0;k<numparams;k++)
		    {
		      // adapt the standard deviation
		      // so that it goes towards 0.3333:
		      rw[k] *= exp((totalacc[k]/100.0-0.3333)*2.0);
		      
		      // initialize the total amount of accpetances again:
		      totalacc[k]=0;
		      
		      if(!silent && k<10)
#ifdef MAIN
			cout << "rw[" << k << "]=" << rw[k] << " ";
#else
			Rcout << "rw[" << k << "]=" << rw[k] << " ";
#endif // MAIN
		    }
		  if(!silent)
#ifdef MAIN
		    cout << endl;
#else
		    Rcout << std::endl;
#endif // MAIN
		}
	      else if(i==(burnin/4))
		{
		  double num_trailing=double(i%100);
		  // traverse the parameters:
		  for(k=0;k<numparams;k++)
		    {
		      // adapt the standard deviation
		      // so that it goes towards 0.3333:
		      rw[k] *= exp((totalacc[k]/num_trailing-0.3333)*2.0);
		      
		      // initialize the total amount of accpetances again:
		      totalacc[k]=0;
		      
		      if(!silent && k<10)
#ifdef MAIN
			cout << "rw[" << k << "]=" << rw[k] << " ";
#else
			Rcout << "rw[" << k << "]=" << rw[k] << " ";
#endif // MAIN
		    }
		  if(!silent)
#ifdef MAIN
		    cout << endl;
#else
		    Rcout << std::endl;
#endif // MAIN
		}
	    }
	}
    }
#ifdef DETAILED_TIMERS
  timers[21][1]=clock();
  timers[21][2]+=(timers[21][1]-timers[21][0]);
#endif // DETAILED_TIMERS

  
#ifdef DETAILED_TIMERS
  timers[22][0]=clock();
#endif // DETAILED_TIMERS

  // **************************************************************
  // MCMC sampling:
  for(i=0;i<numsamples;i++) // traverse the wanted number of samples
    {
      // For each wanted sample, do MCMC iterations 'indep' number
      // of times:
      for(j=0;j<indep;j++)
	{
	  sample_once(data,numrows,numcolumns,
		      numparams, numtemp, params,
		      acc, rw, logprob, T, swaps, tempering_prob,
		      num_hyperparameters, hyper_parameters,
		      log_prior, log_lik, update_prev_log_prob, silent);
	}
      
      if(!silent && i%10==0)
	{
#ifdef MAIN
	  cout << "sample i=" << i+1 << " of " << numsamples << 
	    " j=" << j+1 << " of " << indep << " logprob=" << 
	    logprob[0] << " ";
#else 
	  Rcout << "sample i=" << i+1 << " of " << numsamples << 
	    " j=" <<
	    j+1 << " of " << indep << " logprob=" << 
	    logprob[0] << " ";
#endif // MAIN
	  show_parameter_value(params[0], numparams, param_names, 10);
#ifndef MAIN
	  Rcout << "test" << std::endl;
#endif // MAIN
	}
  
      // copy the result to the return array:
      for(k=0;k<numparams;k++)
	ret[i][k]=params[0][k];
    
    }
#ifdef DETAILED_TIMERS
  timers[22][1]=clock();
  timers[22][2]+=(timers[22][1]-timers[22][0]);
#endif // DETAILED_TIMERS

  // Show tempering swapping info, if wanted:
  if(!silent)
    {
      // Tempering output:
      if(numtemp>1)
	{
#ifdef MAIN
	  printf("Swaps:\n");
	  for(t=0;t<(numtemp-1);t++)
	    printf("%d<->%d: %d\n", t, t+1, swaps[t]);
#else
	  Rcout << "Swaps:" << std::endl;
	  for(t=0;t<(numtemp-1);t++)
	    Rcout << t << "<->" << t+1 << ": " << swaps[t] << std::endl;
#endif // MAIN
	}
    }

  delete [] logprob;
  delete [] acc;
  delete [] totalacc;
  delete [] rw;
  delete [] swaps;
  doubledelete(params,numtemp);

  if(show_graphs)
    {
      for(k=0;k<numparams;k++)
	{
	  double *theta=new double[numsamples];
	  for(i=0;i<numsamples;i++)
	    theta[i]=ret[i][k];
	  
	  show_mcmc_parameter(theta, numsamples, param_names[k],false);
	  delete [] theta;
	}

      for(j=0;j<(numparams-1);j++)
	for(k=j+1;k<numparams;k++)
	  {
	    double *theta1=new double[numsamples];
	    double *theta2=new double[numsamples];
	    for(i=0;i<numsamples;i++)
	      {
		theta1[i]=ret[i][j];
		theta2[i]=ret[i][k];
	      }

#ifdef MAIN
	    show_mcmc_scatter(theta1, theta2, numsamples, 
			      param_names[j],param_names[k]);
#endif // MAIN
	    delete [] theta1;
	    delete [] theta2;
	  }
    }
  
  return ret;
}


// Importance sampling method for calculating model likelihood
// using a multinormal approximation of the posterior samples
// as the proposal distribution.
double log_model_likelihood_multinormal(unsigned int num_importance_samples,
			    
			    double **mcmc_params, 
			    unsigned int numparams, unsigned int numsamples,
			    
			    double **data, unsigned int numrows, unsigned int numcolumns,
			    
			    unsigned int num_hyperparameters, double *hyper_parameters,
			    
			    double (*log_prior)(unsigned int /*numparams*/, 
						double* /* params */,
						unsigned int /* num_hyperparameters */,
						double* /* hyper_parameters */),
			    
			    double (*log_lik)(double** /*data*/, 
					      unsigned int /* numrows */, 
					      unsigned int /* numcolumns */,
					      unsigned int /*numparams*/, 
					      double* /* params */),
			    
			    bool silent,
                            bool dont_keep_cholesky)
{
#ifdef DETAILED_TIMERS
  timers[23][0]=clock();
#endif // DETAILED_TIMERS
  
  int i, num_imp=num_importance_samples;
  
  //Rcout << "log_model_likelihood_multinormal" << std::endl;

  // **************************************************************
  // Importance sampling for finding the marginal data probability
  // Get mean and correlation for the coefficients 
  // in the different models:
      
  // Fetch the moments:
  double *mu_coefs=mean_of_vectors(numsamples, mcmc_params, numparams);
  double **sigma_coefs=estimated_variance(numsamples, mcmc_params,
					  mu_coefs, numparams);
  
  // Initialize GSL random generator:
  //gsl_rng *rptr=gsl_rng_alloc(gsl_rng_rand48);
  //gsl_rng_set(rptr, rand()); 
  
  double *lp_lw=new double[num_imp];
  
  // Importance sampling in order to calculate
  // the model probabilities:
  for(i=0;i<num_imp;i++)
    {
      // sample from a multinormal proposal distribution having 
      // the same moments as the posterior MCMC samples:
      double **csample=multinormal_sample(1, mu_coefs, 
					  sigma_coefs, 
					  numparams, 
					  dont_keep_cholesky);
      if(!csample)
	{
	  delete [] lp_lw;
	  delete [] mu_coefs;
	  doubledelete(sigma_coefs,numparams);
	  //gsl_rng_free(rptr);
#ifdef DETAILED_TIMERS
	  timers[23][1]=clock();
	  timers[23][2]+=(timers[23][1]-timers[23][0]);
#endif // DETAILED_TIMERS
  
	  return MISSING_VALUE;
	}


      // make a proposal parameter structure:
      double *par=*csample;
      
      // calculate the proposal distribution density:
      double lw=pdf_multinormal(par, mu_coefs,
				sigma_coefs, numparams, 
				true);
      if(lw==MISSING_VALUE)
	{
	  delete [] lp_lw;
	  delete [] mu_coefs;
	  doubledelete(sigma_coefs,numparams);
	  //gsl_rng_free(rptr);
#ifdef DETAILED_TIMERS
	  timers[23][1]=clock();
	  timers[23][2]+=(timers[23][1]-timers[23][0]);
#endif // DETAILED_TIMERS
	  return MISSING_VALUE;
	}

      // fetch the logprob for the proposal:
      double lp=log_prob(data,numrows,numcolumns,numparams,par,
			 num_hyperparameters, hyper_parameters,
			 log_prior, log_lik);

      if(lp>-1e+200 && lp<1e+200 && lw>-1e+200 && lw<1e+200)
	lp_lw[i]=lp-lw;
      else
	lp_lw[i]=-1e+200;

      // cleanup:
      doubledelete(csample,1);
    }

  double lp_lw_max=find_statistics(lp_lw,num_imp,MAX);
  
  // probaiblity sum, contribution and prior sum
  double probsum=0.0;
  for(i=0;i<num_imp;i++)
    probsum+=exp(lp_lw[i]-lp_lw_max);

  // Monte Carlo estimate for the marginal data density
  probsum/=((double) num_imp);

  double lprobsum=log(probsum)+lp_lw_max;
  // Show the model marginal data density:
  if(!silent)
    {
#ifdef MAIN
      printf("probsum_0=%g\n", probsum);
      printf("lprobsum=%g\n", lprobsum);
      printf("Id stategy:%d\n",(int) id_strategy);
#else
      Rcout << "probsum_0=" << probsum << std::endl;
      Rcout << "lprobsum =" << lprobsum << std::endl;
      Rcout << "Id strategy =" << id_strategy << std::endl;
#endif // MAIN
    }
  
  delete [] lp_lw;
  delete [] mu_coefs;
  doubledelete(sigma_coefs,numparams);

#ifdef DETAILED_TIMERS
  timers[23][1]=clock();
  timers[23][2]+=(timers[23][1]-timers[23][0]);
#endif // DETAILED_TIMERS
  
  return lprobsum;
}

// show_parameter_value : prints current parameter array to screen:
void show_parameter_value(double *param, unsigned int numparam, char **parnames,
			  int max_params)
{
  unsigned int showmax=numparam;
  if(max_params!=(int)MISSING_VALUE && max_params>0 &&
     numparam>(unsigned int) max_params)
    showmax=max_params;
#ifdef MAIN
  for(unsigned int i=0;i<showmax;i++)
    if(parnames)
      cout << parnames[i] << "=" << param[i] << " ";
    else
      cout << "par" << i+1 << "=" << param[i] << " ";
  cout << endl;
#else
  for(unsigned int i=0;i<showmax;i++)
    if(parnames)
      Rcout << parnames[i] << "=" << param[i] << " ";
    else
      Rcout << "par" << i+1 << "=" << param[i] << " ";
  Rcout << std::endl;
#endif // MAIN
}


// show_parameters: Uses the external programs
// 'vvgraph' and 'histogramme' from the hydrasub
// package to show the sampling for a given
// parameter. If these are not present, this
// methods should be set on silent (only spacing
// between independent samples are shown).
// par: parameter samples
// N: number of samples
// parname: name of parameter
// silent: if set, turns on silent modus
void show_mcmc_parameter(double *par, unsigned int N, char *parname, 
			 bool silent, 
			 char *filestart)
{
  int len=N; 
  // calculate one-step autocorrelation:
  double rho=get_auto_correlation(par, len);
  // number of independent samples, using 
  // AR(1) model and Gelman's formula for independent samples:
  // double n_indep=len/2.0/(0.5+rho/(1.0-rho));
  // spacing=N/n_indep
  double spacing=2.0*(0.5+rho/(1.0-rho));

#ifdef MAIN
  cout << parname << " mean=" << find_statistics(par,N,MEAN);
  cout << " median=" << find_statistics(par,N,MEDIAN);
  cout << " 95% post.cred.int=(" << 
    find_statistics(par,N,PERCENTILE_2_5) << " , " <<
    find_statistics(par,N,PERCENTILE_97_5) << ")" << endl;
#else
  Rcout << parname << " mean=" << find_statistics(par,N,MEAN);
  Rcout << " median=" << find_statistics(par,N,MEDIAN);
  Rcout << " 95% post.cred.int=(" << 
    find_statistics(par,N,PERCENTILE_2_5) << " , " <<
    find_statistics(par,N,PERCENTILE_97_5) << ")" << std::endl;
#endif // MAIN

  // if silent, this is all that is to be done:
  if(silent)
    return;
  
  // Show spacing:
#ifdef MAIN
  cout << parname << " - spacing between independent samples:" << 
    spacing << endl;
#else
  Rcout << parname << " - spacing between independent samples:" << 
    spacing << std::endl;
#endif // MAIN

#ifdef MAIN
  char cmd[1000], filename[1000]; // command string + file name string
  FILE *p; // file pointer
  int i;
  
  // show sampling history using 'vvgraph':
  if(!filestart || !*filestart)
    {
      snprintf(cmd, 999, "vvgraph -x \"%s\" 2> /dev/null > /dev/null", parname);
      p=popen(cmd,"w");
    }
  else
    {
      snprintf(filename, 999, "%s_%s_ts.txt", filestart, parname);
      p=fopen(filename,"w");
    }
  fprintf(p,"# Column 1: %s\n",  parname);
  fprintf(p,"###################\n");
  for(i=0;i<len;i++)
    fprintf(p,"%d %f\n", i+1, par[i]);
  if(!filestart || !*filestart)
    pclose(p);
  else
    fclose(p);

  // show histogram using 'histogramme':
  if(!filestart || !*filestart)
    {
      snprintf(cmd, 999,
	       "histogramme -x \"%s\" -t \"%s\" 2> /dev/null > /dev/null", 
	       parname,  parname); 
      p=popen(cmd, "w");
    }
  else
    {
      snprintf(filename, 999, "%s_%s_hist.txt", filestart, parname);
      p=fopen(filename,"w");
    }
  for(i=0;i<len;i++)
    fprintf(p,"%f\n", par[i]);
  if(!filestart || !*filestart)
    pclose(p);
  else
    fclose(p);
#endif // MAIN
}


#ifdef MAIN
// show_scatter: Shows scatterplots of the samples for
// two parameter.
// par1: parameter samples for parameter 1
// par2: parameter samples for parameter 2
// N: number of samples
// parname1: parameter name 1
// parname2: parameter name 2
void show_mcmc_scatter(double *par1, double *par2, int N, 
		       char *parname1, char *parname2, char *filestart)
{
  FILE *p; // file pointer
  char cmd[1000], filename[1000]; // command string + file name string
  int i, len=N;
  
  // show scatterplot using the 'vvgraph' program:

  if(!filestart || !*filestart)
    {
      snprintf(cmd, 999, "vvgraph -x \"%s\" -y \"%s\"", parname1, parname2);
      p=popen(cmd,"w");
    }
  else
    {
      snprintf(filename, 999, "%s_%s_%s_scatter.txt", filestart, parname1,parname2);
      p=fopen(filename,"w");
    }
  fprintf(p, "# Column 1: %s vs %s\n", parname1, parname2);
  fprintf(p, "# Column 1 - type: dot\n");
  fprintf(p, "#####################\n");
  for(i=0;i<len;i++)
    fprintf(p,"%f %f\n", par1[i], par2[i]);
  if(!filestart || !*filestart)
    pclose(p);
  else
    fclose(p);
}
#endif // MAIN



// *****************************************************
// Necessary structures: HydDateTime
// Pre-requires ctime and HydDate.
// *****************************************************



// ********************************************************
//
// MEASUREMENTS
//
// Structure used for representing the measurements
// found in the data file, containing site, time,
// mean logarithmic val, standard deviation
// of the logarithmic val (can have been smoothed),
// and number of vals contain in this mean (log-)val.
//
// ********************************************************

struct series_measurements
{
  unsigned int serie_num; /* belongs to which series? */
  int  site /*site*/; 
  int  n; /* number of val measurement in this measurement line*/;
  HydDateTime dt; // Time in date/time format (only for recent events)
  double tm /* time=-age */, 
    meanval /* mean (logarithmic) value*/, 
    sd /* standard deviation of (logarithmic) value */;
  int index;

  series_measurements(); 
  
  void print(void);
};

series_measurements::series_measurements()
{
  serie_num=site=n=0; 
  dt=NoHydDateTime; 
  tm=meanval=sd=MISSING_VALUE; 
  index=0;
}
  

class measurement_cluster
{
  void cleanup(void);
public:
  HydDateTime dt; // Time in date/time format (only for recent events)
  double tm /* time=-age */;
  unsigned int num_measurements;

  int *serie_num; /* belongs to which series? */
  int *site /*site*/; 
  int  *n /* number of val measurement in this measurement line*/;
  double *meanval /* mean (logarithmic) value*/, 
    *sd /* standard deviation of (logarithmic) value */;
  int *index; // index in the state vector
  
  double **corr_matrix;
  
  // empty measurement cluster:
  measurement_cluster();
  measurement_cluster(double time_, HydDateTime dtime=NoHydDateTime); 
  // Filled measurement cluster:
  measurement_cluster(double time_, unsigned int number_of_measurements, 
		      int *series, int *sites,
		      int *N, 
		      double *mean_val, double *standard_dev,
		      int *index_, HydDateTime dtime=NoHydDateTime);
  measurement_cluster(measurement_cluster *orig);

  ~measurement_cluster();

  void copy(measurement_cluster *orig);
  void add_measurement(series_measurements &m);

  void print(void);
};


// ********************************************************
//
// PRIOR
//
// Structure that contains info enough to form the
// prior distribution.
// ********************************************************

struct prior
{
  int is_log;
  
  double mu_1, mu_2, mu_m, mu_s;
  // prior for expectancy terms
  
  double dt_1, dt_2, ldt_m, ldt_s;
  // prior for drift terms
  
  double s_1, s_2, ls_m, ls_s;
  // prior for noise terms
  
  double beta_1, beta_2, beta_m, beta_s;
  // prior for regression terms for supplementary series
  // prior for regression term from global full external series
  //       when used for the main series.
  
  
  double os_1, os_2, los_m, los_s;
  // prior for observational noise terms
  
  // Prior for initial values:
  double init_1, init_2, init_m, init_s;
  

  // prior for linear trend
  double lin_1, lin_2, lin_m, lin_s;
  
  double mean_val;

  // correlations terms will have a uniform prior from -0.2 to 1.
  
  void set(double mu1, double mu2, double dt1,   double dt2,
	   double s1,  double s2, double lin1, double lin2, 
	   double beta1, double beta2, 
	   double init1, double init2,
	   double obs_sd1=MISSING_VALUE,double obs_sd2=MISSING_VALUE,
	   double meanval=MISSING_VALUE /* used for setting prior for
					   log-transformed observational noise
					*/);
  void copy(prior *orig);
  
  void show(void);

  // constructors:
  prior(double mu1, double mu2, double dt1,   double dt2,
	double s1,  double s2, double lin1, double lin2,  
	double beta1, double beta2, double init1, double init2);
  prior(double mu1, double mu2, double dt1,   double dt2,
	double s1,  double s2, double lin1, double lin2,  
	double beta1, double beta2, double init1, double init2,
	double obs_sd1,double obs_sd2, double meanval);
  prior(int islog, double mu1, double mu2, double dt1,   double dt2,
	double s1,  double s2, double lin1, double lin2,  
	double beta1, double beta2, double init1, double init2,
	double obs_sd1,double obs_sd2, double meanval);
  prior(prior *orig);
  prior(char *infile,double meanval);
};


// ********************************************************
//
// SERIES
//
// Structure that contains info enough to handle
// a single timeseries. Includes the series_measurements,
// the prior and switches for the structural model.
//
// ********************************************************

class series
{
public:
  series_measurements *meas, *comp_meas;
  unsigned int meas_len, comp_len;
  prior *pr;
  char name[100];
  double mean_time, mean_val; 
  bool linear_time_dep;
  bool no_layers;
  
  // Current parameters used in loglik:
  // (Full model parameters)
  double *mu, **init_mu;
  double **pull;
  double **sigma;
  double *corr;
  double ***paircorr;
  double beta, *lin_t;
  double obs_sd;
  
  // Trigonometric functions:
  unsigned int num_per;
  double Tper[100]; 
  double beta_sin[100], beta_cos[100];
  
  unsigned int numlayers;
  
  int init_treatment, allow_positive_pulls, regional_init, layered_init;
  double init_time, init_value;
  HydDateTime init_dt;
  
  int regional_mu, regional_lin_t,regional_pull[100], regional_sigma[100];
  int no_sigma[100], no_pull_lower;
  int sigma_correlated[100], sigma_pairwise_correlated[100],sigma_1dim[100];
  int indicator_corr[100], indicator_corr2[100], indicator_sigma[100],
    indicator_pull[100], indicator_mu, indicator_lin_t;
  int time_integral[100];
  int serie_num;

#ifdef MAIN
  void read_series(char **&argv, int &argc, bool is_main_series, 
		   int serie_number,
		   bool *use_site,
		   double start_time=MISSING_VALUE, 
		   double end_time=MISSING_VALUE);
#endif // MAIN
  
  void set_series(char *seriename, int serienum, 
		  double *X_time, double *X_value,
		  unsigned int n, unsigned int num_layers, 
		  prior *useprior=NULL, bool is_datetime=false,
		  double *sd=NULL, int *num_meas_per_value=NULL, int *sites=NULL);
  
  void set_series(double **X_time, unsigned int n, unsigned int num_layers);
  
  void cleanup(void);
  void init(void);

  series(); // constructor
  ~series(); // destructor, ready series for deletion
};


// ********************************************************
//
// PARAMS 
//
// Should contained transformed (usually log-transformed) 
// parameter values:
// ********************************************************

class params
{
public:
  int numparam; // number of parameters
  double *param, // transformed parameter value vector
    log_prob, // logarithmic probability (prior times likelihood)
    log_lik; // logarithmic likelihood
  bool iscomplex;
  double cycles[10];

  void cleanup(void); // cleanup the structure
  void copy(params *orig); // copy the contents of another 
  //'params' structure

  // constructors:
  params(); // creates an empty parameter value vector
  params(int numparams); // creates the parameter value vector
  params(double *values, int numparams); // creates the
  // parameter value vector and fills it with values
  params(params *orig); // copies the contents of another 
  //'params' structure

  // destructor:
  ~params();
};




// *******************
// *******************
// Functions:
// *******************
// *******************

// logarithmic factorial:
double lfactorial(int x);

// Draw pairwise correlations as a vector from the prior distribution
double *draw_correlations(int num_sites);

// Make a diffusion/covariance matrix from a correlation vector
double **make_sigma2(double *cor, int num_sites);

// Estimate (by Monte Carlo method) the probability of a correlation drawn 
// from the prior resulting in a positively definite covariance matrix.
double p_positive_correlations(unsigned int num_sites,unsigned int N);


#ifdef MAIN
// Print usage and exit:
void usage(void);
#endif // MAIN


// get_measurements: read the measurement file and store it in
// an array of type 'measurement'. Allows for removing some 
// measurements, by site or by age.
series_measurements *get_measurements(char *infile,
				      unsigned int *len, bool *use_site,
				      double start_time=MISSING_VALUE, 
				      double end_time=MISSING_VALUE);

// get_correlations: read the series correlation file and store it in
// a timeseries of arrays, corresponding to corr_1_2, corr_1_3, ... ,
// coor_1_S, corr_2_3, ..., corr_2_S, ...,  corr_S-1_S. In total S*(S-1)/2
// correlations, where S is the number of series. Allows for removing some 
// measurements, by site or by age.
double **get_correlations(char *infile, int *len, int **sites, 
			  double **t, bool datetime_format,
			  bool *use_site,
			  double start_time=MISSING_VALUE, 
			  double end_time=MISSING_VALUE);


// transform_paramete: Transforming from original 
// parametrization to transformed parameters.
double transform_parameter(double val, transform_type type);

// invtransform_paramete: Transforming from transformed
// parametrization to original parameters.
double invtransform_parameter(double val, transform_type type);



// loglik: Returns log-likelihood unless no parameteers
// are given or return_correlation is switched on.
// If no parameter array is given (pars=NULL), then it returns 
// the number of parameters. If return_correlation!=0, then
// the correlation function value for the current parameter
// array, site 'return_site' and the time difference 
// 'return_t' is returned. return_correlation, return_site 
// and return_t has default values, see the header file.
double loglik_complex(double *pars, int dosmooth=0, int do_realize=0, 
		      int residual_analysis=0, char *res_filestart=NULL, 
		      int debug=0, char **simulation_files=NULL);
double loglik_real(double *pars, int dosmooth=0, int do_realize=0, 
		   int residual_analysis=0, char *res_filestart=NULL, 
		   int debug=0, char **simulation_files=NULL);
double loglik(double *pars, int dosmooth=0, int do_realize=0, 
	      int residual_analysis=0, char *res_filestart=NULL, 
	      int debug=0, char **simulation_files=NULL,
	      int return_residuals=0,
	      double **residuals_time=NULL, double ***residuals=NULL,
	      double ***prior_expected_values=NULL,
	      int *resid_numcolumns=NULL, int *resid_len=NULL);
// For ML purposes:
double minusloglik(double *pars);

// Analyze residuals:
void analyze_residuals(double *residuals, double *tm, HydDateTime *dt, int len,
		       char *filestart);

// **************************************************************
// logprob: Returns log-probability (prior times likelihood).
// Input: transformed parameter vector, measurements, prior, 
// and 'temperature' (used for tempering purposes).
// **************************************************************
double logprob(params &par, double T, int dosmooth=0, 
	       int do_realize=0, int doprint=0);

// newsample: Sample one new parameter set, using 
// Random Walk Metropolis and parallell tempering on
// the previous sample.
// sample: starts as the previous sample, ends as the
//         current sample. It's an array, since you can
//         have several tempering chains.
// logprob: start as previous logprobs 
//          (log(prior*likelihood) for each chain)
//          and ends as the current logprobs.
// acc: acceptance indicator array. Is set to one for 
//      each parameter that has it's new value accepted.
//      Used for the adaptive phase of the burn-in.
// rw: Random Walk standard deviations for the proposal
//     density of the Metropolis algorithm.
// T: temperature array
// numtemp: Number of tempering chains
// swaps: Incremented each time a swapping of tempering
//        chains takes place.
void newsample(params *sample, double *log_prob, params *acc, params *rw, 
	       double *T, unsigned int numtemp, int *swaps, int dosmooth=0,
	       int do_realize=0);


// mcmc: Perform Markov Chain Monte Carlo sampling.
// numsamples: number of samples.
// burnin: burn-in period.
// indep: number of MCMC iterations between each
//        sample fetched.
// numtemp: Number of tempering chains.
params *layer_mcmc(unsigned int numsamples, unsigned int burnin,
		   unsigned int indep, unsigned int numtemp=1,
		   bool do_importance=true,
		   int dosmooth=0, int do_realization=0,  
		   char *realization_file_start=NULL, double ****x_=NULL, 
		   double *startpar=NULL, double T_ground=1.5,
		   double *model_loglik=NULL, 
		   double *model_dic1=NULL,
		   double *eff_num_param1=NULL, 
		   double *model_dic2=NULL,
		   double *eff_num_param2=NULL);


// show_parameters: Uses the external programs
// 'vvgraph' and 'histogramme' from the hydrasub
// package to show the sampling for a given
// parameter. If these are not present, this
// methods should be set on silent (only spacing
// between independent samples are shown).
// par: parameter samples
// N: number of samples
// parname: name of parameter
// silent: if set, turns on silent modus
void show_parameter(double *par, int N, char *parname, 
		    int silent=0, char *filestart=NULL);

#ifdef MAIN
// show_scatter: Shows scatterplots of the samples for
// two parameter.
// par1: parameter samples for parameter 1
// par2: parameter samples for parameter 2
// N: number of samples
// parname1: parameter name 1
// parname2: parameter name 2
void show_scatter(double *par1, double *par2, int N, 
		  char *parname1, char *parname2, char *filestart=NULL);
#endif // MAIN






#ifdef __linux__  
#include <flexiblas/flexiblas_api.h>
#endif // __linux__

// *************************************************
// Purging mechanism for global variables, so that
// they can be reused.
// *************************************************

bool detailed_loglik=false;

void reset_global_variables(void)
{
  detailed_loglik=false;
  
#ifdef DETAILED_TIMERS
  for(int i=0;i<100;i++)
    for(int j=0;j<3;j++)
      timers[i][j]=0;
#endif // DETAILED_TIMERS


#ifdef __linux__  
  mkl_set_num_threads(1);
#endif // __linux__
  
  
  if(par_name)
    doubledelete(par_name,LARGE_ENOUGH_ARRAY);
  par_name=NULL;
  
  if(par_type)
    delete [] par_type;
  par_type=NULL;
  
  if(par_region)
    delete [] par_region;
  par_region=NULL;
  
  if(par_layer)
    delete [] par_layer;
  par_layer=NULL;
  
  if(par_series)
    delete [] par_series;
  par_series=NULL;

  if(par_trans_type)
    delete [] par_trans_type;
  par_trans_type=NULL;

  if(ser)
    delete [] ser;
  ser=NULL;

  if(corr_from_series)
    delete [] corr_from_series;
  corr_from_series=NULL;

  if(corr_to_series)
    delete [] corr_to_series;
  corr_to_series=NULL;

  if(corr_from_layer)
    delete [] corr_from_layer;
  corr_from_layer=NULL;

  if(corr_to_layer)
    delete [] corr_to_layer;
  corr_to_layer=NULL;

  if(corr_from_index)
    delete [] corr_from_index;
  corr_from_index=NULL;

  if(corr_to_index)
    delete [] corr_to_index;
  corr_to_index=NULL;

  if(series_corr)
    delete series_corr;
  series_corr=NULL;

  if(indicator_array)
    delete [] indicator_array;
  indicator_array=NULL;

  if(meas_tot)
    delete [] meas_tot;
  meas_tot=NULL;

  if(meas_smooth)
    delete [] meas_smooth;
  meas_smooth=NULL;
  
  if(obs_corr_t)
    delete [] obs_corr_t;
  obs_corr_t=NULL;
  
  if(obs_corr)
    doubledelete(obs_corr, obs_corr_len);
  obs_corr=NULL;
  
  if(obs_corr_site)
    delete [] obs_corr_site;
  obs_corr_site=NULL;

  if(x_k_s_kept)
    doubledelete(x_k_s_kept, len_x_k_s_kept);
  x_k_s_kept=NULL;
  
  if(P_k_s_kept)
    tripledelete(P_k_s_kept,len_x_k_s_kept,size_x_k_s_kept);
  P_k_s_kept=NULL;

  if(t_k_smooth)
    delete [] t_k_smooth;
  t_k_smooth=NULL;

  if(dt_k_smooth)
    delete [] dt_k_smooth;
  dt_k_smooth=NULL;

  if(x_k_realized)
    doubledelete(x_k_realized,meas_smooth_len);
  x_k_realized=NULL;
  
  if(x_k_realized_all && numit_realizations_made>0)
    tripledelete(x_k_realized_all,numit_realizations_made, meas_smooth_len);
  x_k_realized_all=NULL;
  numit_realizations_made=0;
  
  if(feed_from_series)
    delete [] feed_from_series;
  feed_from_series=NULL;
  
  if(feed_to_series)
    delete [] feed_to_series;
  feed_to_series=NULL;

  if(feed_from_layer)
    delete [] feed_from_layer;
  feed_from_layer=NULL;

  if(feed_to_layer)
    delete [] feed_to_layer;
  feed_to_layer=NULL;

  if(feed_symmetric)
    delete [] feed_symmetric;
  feed_symmetric=NULL;

  if(beta_feed)
    delete [] beta_feed;
  beta_feed=NULL;
  
  if(extdata)
    delete [] extdata;
  extdata=NULL;

  len_x_k_s_kept=0;
  size_x_k_s_kept=0;
  ext_len=0;
  is_complex=false;
  num_series_feed=0;
  nodata=false;
  numit_realization=100;
  numpar=0;
  num_states=1;
  obs_corr_len=0;
  meas_tot_len=0;
  meas_smooth_len=0;
  useext=0;
  use_indicator=0;
  num_smooth=10; 
  num_series=0;
  real_strat=NO_CENSORING;
  id_strategy=ID_NONE;
  numsites=0;
  num_tot_layers=0;
  some_pairwise_correlations=false;
  p_pos_site_sigma2=0.0;
  p_pos_series_sigma2=0.0;
  num_series_corr=0;
}





// logarithmic factorial:
double lfactorial(int x)
{
  return lgamma(double(x+1));
}

// Make a diffusion/covariance matrix from a correlation vector
double **make_site_sigma2(double *cor, int num_sites)
{
  double **sigma2=Make_matrix(num_sites,num_sites);
  int i,j,p=num_sites;
  
  for(i=0;i<p;i++)
    sigma2[i][i]=1.0;
  
  int k=0;
  for(i=0;i<p;i++)
    for(j=(i+1);j<p;j++)
      sigma2[i][j]=sigma2[j][i]=cor[k++];
  
  return sigma2;
}


double logprior_site_corr(unsigned int numparams, 
			  double *params, // logit-transformed correlations
			  unsigned int /* num_hyperparameters */,
			  double * /*hyper_parameters*/ )
{
  double lp=0.0;
  for(unsigned int i=0;i<numparams;i++)
    lp+= (-0.5*log(2.0*M_PI)-log(2.0)-0.5*params[i]*params[i]/4.0);
  return lp;
}

double loglik_site_corr(double ** /* data */, 
			unsigned int numrows, // stand in for num_sites 
			unsigned int /* numcolumns */,
			unsigned int numparams, 
			double *params)
{
  unsigned int i,j;
  unsigned int num_sites=numrows;
  double *cor=new double[numparams];
  
  for(i=0;i<numparams;i++)
    cor[i]=-1.0+2.0*exp(params[i])/(1.0+exp(params[i]));
  
  double **sigma2=make_site_sigma2(cor,num_sites);
  //Rcout << "loglik_site_corr" << std::endl;
  
  double *e=double_eigenvalues(sigma2, num_sites);
  int isok=1;
  for(j=0;j<num_sites && isok;j++)
    if(e[j]<=0.0)
      isok=0;
  
  delete [] cor;
  doubledelete(sigma2,num_sites);
  delete [] e;

  if(isok)
    return(0.0);
  else
    return(-1e+200);
}		    

void init_site_corr(unsigned int numparams,
		    double *params, 
		    unsigned int /* num_hyperparameters */,
		    double* /* hyperparameters */)
{
  for(unsigned int i=0;i<numparams;i++)
    params[i]=0.0;
}


// Estimate (by MCMC+importance sampling) the probability of a 
// correlation drawn from the prior resulting in a positively 
// definite covariance matrix for the sites in any given combination
// of series and layer.
double p_positive_sitecorr(unsigned int num_sites,unsigned int N)
{
  unsigned int dim=num_sites*(num_sites-1)/2;
  double T=1.0;
  double **params=general_mcmc(N,N,9,1,NULL,num_sites,1,dim,NULL,
			       &T, 0.1, 0, NULL,
			       init_site_corr, logprior_site_corr, 
			       loglik_site_corr, true);

  double lbml=
    log_model_likelihood_multinormal(10*N, params, dim, N,
				     NULL, num_sites,1,
				     0,NULL,
				     logprior_site_corr,loglik_site_corr, 
				     true);
  
#ifdef MAIN
  cout << "lmbl_site=" << lbml << endl;
#else
  Rcout << "lmbl_site=" << lbml << std::endl;
#endif //U MAIN
  

  doubledelete(params,N);

  return exp(lbml);
}

// Make a diffusion/covariance matrix from a correlation vector
double **make_series_sigma2(double *cor)
{
  double **sigma2=Make_matrix(num_tot_layers,num_tot_layers);
  unsigned int i,p=num_series_corr;

  for(i=0;i<num_tot_layers;i++)
    sigma2[i][i]=1.0;
  
  for(i=0;i<p;i++)
    sigma2[corr_from_index[i]][corr_to_index[i]]=
      sigma2[corr_to_index[i]][corr_from_index[i]]=cor[i];

  return sigma2;
}


double logprior_series_corr(unsigned int numparams, 
			    double *params, // logit-transformed correlations
			    unsigned int /* num_hyperparameters */,
			    double * /*hyper_parameters*/ )
{
  double lp=0.0;
  for(unsigned int i=0;i<numparams;i++)
    lp+= (-0.5*log(2.0*M_PI)-log(2.0)-0.5*params[i]*params[i]/4.0);
  return lp;
}

double loglik_series_corr(double ** /* data */, 
			  unsigned int /* numrows */, 
			  unsigned int /* numcolumns */,
			  unsigned int /* numparams */, 
			  double *params)
{
  unsigned int i,j;
  double *cor=new double[num_series_corr];
  
  for(i=0;i<num_series_corr;i++)
    cor[i]=-1.0+2.0*exp(params[i])/(1.0+exp(params[i]));
  
  double **sigma2=make_series_sigma2(cor);
  //Rcout << "loglik_series_corr" << std::endl;
  double *e=double_eigenvalues(sigma2, num_tot_layers);
  int isok=1;
  for(j=0;j<num_tot_layers && isok;j++)
    if(e[j]<=0.0)
      isok=0;
  
  delete [] cor;
  doubledelete(sigma2,num_tot_layers);
  delete [] e;

  if(isok)
    return(0.0);
  else
    return(-1e+200);
}		    

void init_series_corr(unsigned int numparams,
		      double *params, 
		      unsigned int /* num_hyperparameters */,
		      double* /* hyperparameters */)
{
  for(unsigned int i=0;i<num_series_corr;i++)
    params[i]=0.0;
}


// Estimate (by MCMC+importance sampling) the probability of a 
// correlation drawn from the prior resulting in a positively 
// definite covariance matrix for the sites in any given combination
// of series and layer.
double p_positive_seriescorr(unsigned int N)
{
  double T=1.0;
  double **params=general_mcmc(N,N,9,1,
			       NULL,0,0,
			       num_series_corr,NULL,
			       &T, 0.1, 0, NULL, init_series_corr,
			       logprior_series_corr, 
			       loglik_series_corr, true);

  double lbml=
    log_model_likelihood_multinormal(10*N, 
				     params, num_series_corr, N,
				     NULL, 0,0,
				     0,NULL,
				     logprior_series_corr,
				     loglik_series_corr, 
				     true);
  
#ifdef MAIN
  cout << "lmbl_series=" << lbml << endl;
#else
  Rcout << "lmbl_series=" << lbml << std::endl;
#endif // MAIN

  doubledelete(params,N);
  
  return exp(lbml);
}



// ********************************************************
//
// PARAMS 
//
// Should contained transformed (usually log-transformed) 
// parameter values:
// ********************************************************

// Constructors 
params::params()
{
  numparam=0;
  param=NULL;
}

params::params(int numparams)
{
  numparam=numparams;
  param=new double[numparams];
  for(int i=0;i<numparam;i++)
    param[i]=0.0;
}

params::params(params *orig)
{
  numparam=0;
  param=NULL;
  copy(orig);
}

void params::copy(params *orig)
{
  cleanup();
  
  int i;
  numparam=orig->numparam;
  if(orig->param)
    {
      param=new double[numparam];
      for(i=0;i<numparam;i++)
	param[i]=orig->param[i];
    }
}

params::params(double *values, int numparams)
{
  numparam=numparams;
  param=new double[numparam];
  for(int i=0;i<numparam;i++)
    param[i]=values[i];
}

void params::cleanup(void)
{
  if(param)
    delete [] param;
  param=NULL;
  numparam=0;
}


// destructor:
params::~params()
{
  cleanup();
}



// ********************************************************
//
// PRIOR
//
// Structure that contains info enough to form the
// prior distribution.
// ********************************************************

// constructors:
prior::prior(double mu1, double mu2, double dt1, double dt2,
	     double s1, double s2, double lin1, double lin2, 
	     double beta1, double beta2, double init1, double init2)
{
  is_log=0;
  set(mu1,mu2,dt1,dt2,s1,s2,lin1,lin2,beta1,beta2, init1, init2);
}

prior::prior(double mu1, double mu2, double dt1, double dt2,
	     double s1, double s2, double lin1, double lin2, 
	     double beta1, double beta2, double init1, double init2,
	     double obs_sd1,double obs_sd2,double meanval)
{
  is_log=0;
  set(mu1,mu2,dt1,dt2,s1,s2,lin1,lin2,beta1,beta2, 
      init1, init2,obs_sd1,obs_sd2,meanval);
}

prior::prior(int islog, double mu1, double mu2, double dt1, double dt2,
	     double s1, double s2, double lin1, double lin2, 
	     double beta1, double beta2, double init1, double init2,
	     double obs_sd1,double obs_sd2,double meanval)
{
  is_log=islog;
  set(mu1,mu2,dt1,dt2,s1,s2,lin1,lin2,beta1,beta2, 
      init1, init2,obs_sd1,obs_sd2,meanval);
}

prior::prior(prior *orig)
{
  copy(orig);
}

void prior::copy(prior *orig)
{
  if(!orig)
    return;

  is_log=orig->is_log;
  set(orig->mu_1,orig->mu_2,orig->dt_1,orig->dt_2,orig->s_1,orig->s_2,
      orig->lin_1,orig->lin_2,orig->beta_1,orig->beta_2,
      orig->init_1, orig->init_2, orig->os_1,orig->os_2,
      orig->mean_val);
}

void prior::show(void)
{
  char str[1000];
  snprintf(str,999,"mu  in (%7.3f,%7.3f),  m =%7.3f s =%7.3f),\n"
	  "dt  in (%7.3f,%7.3f),  lm=%7.3f ls=%7.3f),\n"
	  "s   in (%7.3f,%7.3f),  lm=%7.3f ls=%7.3f),\n"
	  "lin in (%7.3f,%7.3f),  m =%7.3f s =%7.3f),\n" 
	  "beta in (%7.3f,%7.3f), m =%7.3f s =%7.3f),\n"
	  "init in (%7.3f,%7.3f), m =%7.3f s =%7.3f),\n"
	  "obs in (%7.3f,%7.3f),  lm=%7.3f ls=%7.3f).\n"
	  "is_log=%d.",
	  mu_1,mu_2,mu_m,mu_s,
	  dt_1,dt_2,ldt_m,ldt_s,
	  s_1,s_2,ls_m,ls_s,
	  lin_1,lin_2,lin_m,lin_s,
	  beta_1,beta_2,beta_m,beta_s,
	  init_1,init_2,init_m,init_s,
	  os_1,os_2,los_m,los_s,
	  is_log);

#ifdef MAIN
  cout << str << endl;
#else
  Rcout << str << std::endl;
#endif // MAIN
}

void prior::set(double mu1, double mu2, double dt1, double dt2,
		double s1, double s2, double lin1, double lin2, 
		double beta1, double beta2, double init1, double init2,
		double obs_sd1,double obs_sd2,double meanval)
{
  mu_1=mu1;
  mu_2=mu2;
  dt_1=dt1;
  dt_2=dt2;
  s_1=s1;
  s_2=s2;
  lin_1=lin1;
  lin_2=lin2;
  beta_1=beta1;
  beta_2=beta2;
  init_1=init1;
  init_2=init2;
  mean_val=meanval;
  
  mu_m=(mu_1+mu_2)/2.0;
  mu_s=(mu_2-mu_1)/2.0/1.96;
  ldt_m=(log(dt_1)+log(dt_2))/2.0;
  ldt_s=(log(dt_2)-log(dt_1))/2.0/1.96;
  ls_m=(log(s_1)+log(s_2))/2.0;
  ls_s=(log(s_2)-log(s_1))/2.0/1.96;
  lin_m=(lin_1+lin_2)/2.0;
  lin_s=(lin_2-lin_1)/2.0/1.96;
  beta_m=(beta_1+beta_2)/2.0;
  beta_s=(beta_2-beta_1)/2.0/1.96;
  init_m=(init_1+init_2)/2.0;
  init_s=(init_2-init_1)/2.0/1.96;
  if(obs_sd1!=MISSING_VALUE && obs_sd2!=MISSING_VALUE)
    {
      if(is_log!=1)
	{
	  os_1=obs_sd1;
	  os_2=obs_sd2;
	}
      else
	{
	  os_1=sqrt(log(1.0+sqrt(1.0+4.0*obs_sd1*obs_sd1/meanval/meanval))-
		    log(2.0));
	  os_2=sqrt(log(1.0+sqrt(1.0+4.0*obs_sd2*obs_sd2/meanval/meanval))-
		    log(2.0));
	}
      if(!silent)
#ifdef MAIN
	cout << os_1 << " " << os_2 << endl;
#else
        Rcout << os_1 << " " << os_2 << std::endl;
#endif // MAIN
	
      los_m=(log(os_1)+log(os_2))/2.0;
      los_s=(log(os_2)-log(os_1))/2.0/1.96;
    }
  else
    os_1=os_2=los_m=los_s=MISSING_VALUE;
}

prior::prior(char *infile, double meanval)
{
  csv c(infile);
  
  mean_val=meanval;

  if(c.get_length()==13)
    {
      int islog=(int) c.get_value(0);
      double mu1=c.get_value(1);
      double mu2=c.get_value(2);
      double dt1=c.get_value(3);
      double dt2=c.get_value(4);
      double s1=c.get_value(5);
      double s2=c.get_value(6);
      double lin1=c.get_value(7);
      double lin2=c.get_value(8);
      double beta1=c.get_value(9);
      double beta2=c.get_value(10);
      double init1=c.get_value(11);
      double init2=c.get_value(12);

      is_log=islog;
      set(mu1,mu2,dt1,dt2,s1,s2,lin1,lin2,beta1,beta2,init1,init2);
    }
  else if(c.get_length()==15)
    {
      int islog=(int) c.get_value(0);
      double mu1=c.get_value(1);
      double mu2=c.get_value(2);
      double dt1=c.get_value(3);
      double dt2=c.get_value(4);
      double s1=c.get_value(5);
      double s2=c.get_value(6);
      double lin1=c.get_value(7);
      double lin2=c.get_value(8);
      double beta1=c.get_value(9);
      double beta2=c.get_value(10);
      double init1=c.get_value(11);
      double init2=c.get_value(12);
      double obs_sd1=c.get_value(13);
      double obs_sd2=c.get_value(14);

      is_log=islog;
      set(mu1,mu2,dt1,dt2,s1,s2,lin1,lin2,beta1,beta2,init1,init2,
	  obs_sd1,obs_sd2,meanval);
    }
  else
    {
#ifdef MAIN
      cerr << "Unknown format for prior file \"" << infile << "\"!" << endl;
      exit(0);
#else
      Rcout << "Unknown format for prior file \"" << infile << "\"!" << std::endl;
#endif // MAIN
    }
}



// Compare measurments, for sorting on time
int compare_meas(const void *i, const void *j) 
{
  double v =  ((series_measurements *)i)->tm - ((series_measurements *)j)->tm;
  if(v < 0)
    return -1;
  else if(v == 0)
    return 0;
  else
    return 1;
} /* compar */

// Compare measurments, for sorting on time
int compare_meas_cluster(const void *i, const void *j) 
{
  double v =  ((measurement_cluster *)i)->tm - ((measurement_cluster *)j)->tm;
  if(v < 0)
    return -1;
  else if(v == 0)
    return 0;
  else
    return 1;
} /* compar */

void series_measurements::print(void)
{
#ifdef MAIN  
  if(sd!=MISSING_VALUE && n!=(int) MISSING_VALUE)
    printf("ser=%d site=%d tm=%9.5f val=%9.5f sd=%9.5f n=%d",
	   serie_num,site,tm,meanval,sd,n);
  else if(sd!=MISSING_VALUE)
    printf("ser=%d site=%d tm=%9.5f val=%9.5f sd=%9.5f",
    	   serie_num,site,tm,meanval,sd);
  else
    printf("ser=%d site=%d tm=%9.5f val=%9.5f",
    	   serie_num,site,tm,meanval);
  if(dt!=NoHydDateTime)
    printf(" dt=%s", dt.syCh(1));
  printf("\n");
#else
  if(sd!=MISSING_VALUE && n!=(int) MISSING_VALUE)
    Rcout << "ser=" << serie_num << " site=" << site << " tm=" << tm <<
      " val=" << meanval << " sd=" << sd << " n=" << n;
  else if(sd!=MISSING_VALUE)
    Rcout << "ser=" << serie_num << " site=" << site << " tm=" << tm <<
      " val=" << meanval << " sd=" << sd;
  else
    Rcout << "ser=" << serie_num << " site=" << site << " tm=" << tm <<
      " val=" << meanval;
  if(dt!=NoHydDateTime)
    Rcout << " dt=" << dt.syCh(1);
  Rcout << std::endl;
#endif //MAIN  
}

// ********************************************************
//
// Measurement_cluster 
//
// Represents a set of measurements stemming from the
// same time point.
// ********************************************************

measurement_cluster::measurement_cluster()
{
  tm=MISSING_VALUE;
  dt=NoHydDateTime;
  num_measurements=0;

  serie_num=site=NULL;
  n=index=NULL;
  meanval=sd=NULL;
  corr_matrix=NULL;
}

measurement_cluster::measurement_cluster(double time_, HydDateTime dtime)
{
  tm=time_;
  dt=dtime;
  num_measurements=0;

  serie_num=site=NULL;
  n=index=NULL;
  meanval=sd=NULL;
  corr_matrix=NULL;
}

// Filled measurement cluster:
measurement_cluster::measurement_cluster(double time_,
					 unsigned int number_of_measurements, 
					 int *series,
					 int *sites, int *N, 
					 double *mean_val, double *standard_dev,
					 int *index_, HydDateTime dtime)
{
  tm=time_;
  dt=dtime;
  num_measurements=number_of_measurements;

  if(num_measurements>0)
    {
      serie_num=new int[num_measurements];
      site=new int[num_measurements];
      n=new int[num_measurements];
      index=new int[num_measurements];
      meanval=new double[num_measurements];
      sd=new double[num_measurements];

      for(unsigned int i=0;i<num_measurements;i++)
	{
	  if(series)
	    serie_num[i]=series[i];
	  if(sites)
	    site[i]=sites[i];
	  if(N)
	    n[i]=N[i];
	  if(index_)
	    index[i]=index_[i];
	  if(meanval)
	    meanval[i]=mean_val[i];
	  if(standard_dev)
	    sd[i]=standard_dev[i];
	}
    }
  else
    {
      serie_num=site=NULL;
      n=index=NULL;
      meanval=sd=NULL;
    }
  corr_matrix=NULL;
}

measurement_cluster::measurement_cluster(measurement_cluster *orig)
{
  tm=MISSING_VALUE;
  dt=NoHydDateTime;
  num_measurements=0;

  serie_num=site=NULL;
  n=index=NULL;
  meanval=sd=NULL;
  corr_matrix=NULL;

  copy(orig);
}

void measurement_cluster::cleanup(void)
{
  if(num_measurements>0)
    {
      if(serie_num)
	delete [] serie_num;
      serie_num=NULL;
      if(site)
	delete [] site;
      site=NULL;
      if(n)
	delete [] n;
      n=NULL;
      if(index)
	delete [] index;
      index=NULL;
      if(meanval)
	delete [] meanval;
      meanval=NULL;
      if(sd)
	delete [] sd;
      sd=NULL;
      if(corr_matrix)
	doubledelete(corr_matrix,num_measurements);
      corr_matrix=NULL;
    }
  num_measurements=0;
  tm=MISSING_VALUE;
}

measurement_cluster::~measurement_cluster()
{
  cleanup();
}

void measurement_cluster::copy(measurement_cluster *orig)
{
  tm=orig->tm;
  dt=orig->dt;
  num_measurements=orig->num_measurements;
  unsigned int i,j;

  if(num_measurements>0)
    {
      serie_num=new int[num_measurements];
      site=new int[num_measurements];
      n=new int[num_measurements];
      index=new int[num_measurements];
      meanval=new double[num_measurements];
      sd=new double[num_measurements];

      for(i=0;i<num_measurements;i++)
	{
	  if(orig->serie_num)
	    serie_num[i]=orig->serie_num[i];
	  if(orig->site)
	    site[i]=orig->site[i];
	  if(orig->n)
	    n[i]=orig->n[i];
	  if(orig->index)
	    index[i]=orig->index[i];
	  if(orig->meanval)
	    meanval[i]=orig->meanval[i];
	  if(orig->sd)
	    sd[i]=orig->sd[i];
	}

      if(orig->corr_matrix)
	{
	  corr_matrix=Make_matrix(num_measurements,num_measurements);
	  for(i=0;i<num_measurements;i++)
	    for(j=0;j<num_measurements;j++)
	      corr_matrix[i][j]=orig->corr_matrix[i][j];
	}
      else
	corr_matrix=NULL;
    }
  else
    {
      serie_num=site=NULL;
      n=index=NULL;
      meanval=sd=NULL;
      corr_matrix=NULL;
    }
}

void measurement_cluster::add_measurement(series_measurements &m)
{
  unsigned int i;

  if(num_measurements>0)
    {
      if(m.tm!=tm)
	{
#ifdef MAIN
	  printf("Time mismatch, when adding measurement:\n");
	  m.print();
	  printf("...to a measurement cluster with time %lf\n", tm);
	  exit(0);
#else
	  Rcout << "Time mismatch, when adding measurement!" << std::endl;
#endif // MAIN
	}

      for(i=0;i<num_measurements;i++)
	if(m.site==site[i] && (int)m.serie_num==serie_num[i])
	  {
#ifdef MAIN
	    printf("Added a new measurement which was already included:\n");
	    m.print();
	    printf("%d measurement cluster with time %lf\n", i, tm);
	    exit(0);
#else
	    Rcout << "Added a new measurement which was already included!"
		  << std::endl;
#endif // MAIN
	  }
    }

  unsigned int num=num_measurements+1;
  int *index_=new int[num];
  int *series=new int[num];
  int *sites=new int[num];
  int *N=new int[num];
  double *mm=new double[num];
  double *ss=new double[num];
  
  for(i=0;i<num_measurements;i++)
    {
      if(index)
	index_[i]=index[i];
      else
	index_[i]=(int) MISSING_VALUE;

      if(serie_num)
	series[i]=serie_num[i];
      else
	series[i]=(int) MISSING_VALUE;

      if(site)
	sites[i]=site[i];
      else
	sites[i]=(int) MISSING_VALUE;

      if(n)
	N[i]=n[i];
      else
	N[i]=(int) MISSING_VALUE;

      if(meanval)
	mm[i]=meanval[i];
      else
	mm[i]=MISSING_VALUE;

      if(sd)
	ss[i]=sd[i];
      else
	ss[i]=MISSING_VALUE;
    }
  
  index_[i]=m.index;
  series[i]=m.serie_num;
  sites[i]=m.site;
  N[i]=m.n;
  mm[i]=m.meanval;
  ss[i]=m.sd;

  cleanup();
  
  num_measurements=num;
  index=index_;
  serie_num=series;
  site=sites;
  n=N;
  meanval=mm;
  sd=ss;

  tm=m.tm;
  dt=m.dt;
}


void measurement_cluster::print(void)
{
#ifdef MAIN
  printf("tm=%9.5f #measurements=%d\n", tm, num_measurements);
  for(unsigned int i=0;i<num_measurements;i++)
    {
      if(sd && sd[i]!=MISSING_VALUE && n && n[i]!=(int) MISSING_VALUE)
	printf("ser=%d site=%d val=%9.5f sd=%9.5f n=%d\n",
	       serie_num[i],site[i],meanval[i],sd[i],n[i]);
      else if(sd && sd[i]!=MISSING_VALUE)
	printf("ser=%d site=%d val=%9.5f sd=%9.5f\n",
	       serie_num[i],site[i],meanval[i],sd[i]);
      else
	printf("ser=%d site=%d val=%9.5f\n",
	       serie_num[i],site[i],meanval[i]);
    }
  if(corr_matrix)
    {
      printf("Correlation matrix:\n");
      for(unsigned int i=0;i<num_measurements;i++)
	{
	  for(unsigned int j=0;j<num_measurements;j++)
	    printf("%7.5lf ", corr_matrix[i][j]);
	  printf("\n");
	}
    }
  if(dt!=NoHydDateTime)
    printf(" dt=%s", dt.syCh(1));
  printf("\n");
#endif // MAIN
}


// ********************************************************
// get_measurements: read the measurement file and store 
// it in an array of type 'measurement'. Allows for 
// removing some measurements, by site or by age.
// ********************************************************

bool nosites=false;
series_measurements *get_measurements(char *infile,
				      unsigned int *len, 
				      bool *use_site,
				      double start_time, double end_time)
{
  // buffer variables for the contents of each line:
  double tm,meanval,sd; 
  int site,n;
  int yr,mnt,day,hour,min,sec; // only for recent time series
  // (time will be counted in minutes rather than m,illions of years)
  static int first=1;
  
  unsigned int i,lines=0; // index variables
  
  ifstream in;
  char str[10000];

  in.open(infile, ios::in);
  if(in.fail())
    {
#ifdef MAIN
      printf("Couldn't open file \"%s\"", infile);
      exit(0);
#else
      Rcout << "Couldn't open file \"" << infile << "\"!" << std::endl;
#endif // MAIN
    }
  
  in.getline(str,9999);
  while(!in.eof())
    {
      // attempt to read the contents:
      if(sscanf(str,"%04d-%02d-%02d %02d:%02d:%02d;%lf", 
		     &yr,&mnt,&day,&hour,&min,&sec, 
		     &meanval)==7)
	{
	  HydDateTime dt(yr,mnt,day,hour,min,sec),start=dt.StartOfYear(),ref(1970);
	  int secs=dt-start;
	  double tm=(double)yr + double(secs)/
	    double(dt.daysInYear())/1440.0/60.0;
	  if((start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    // increase the number of lines to be read
	    lines++;
	}
      else if(!nosites && sscanf(str,"%d %lf %lf %lf %d", &site, 
                &tm, &meanval,&sd,&n)==5)
	{
	  // if this went well and the measurement is not to be removed...
	  if(use_site[site] && (start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    // increase the number of lines to be read
	    lines++;
	}
      else if(!nosites && sscanf(str,"%d %lf %lf %lf", &site, 
		     &tm, &meanval,&sd)==4)
	{
	  if(use_site[site] && (start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    // increase the number of lines to be read
	    lines++;
	}
      else if(!nosites && sscanf(str,"%d %lf %lf", &site,&tm, &meanval)==3 &&
	      getnextint(str,&site)[0]==' ')
	{
	  if(use_site[site] && (start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    // increase the number of lines to be read
	    lines++;
	}
      else if(sscanf(str,"%lf %lf %lf %n", &tm, &meanval, &sd, &n)==4)
	{
	  if((start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    // increase the number of lines to be read
	    lines++;
	}
      else if(sscanf(str,"%lf %lf %lf", &tm, &meanval, &sd)==3)
	{
	  if((start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    // increase the number of lines to be read
	    lines++;
	}
      else if(sscanf(str,"%lf %lf", &tm, &meanval)==2)
	{
	  if((start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    // increase the number of lines to be read
	    lines++;
	}
      else if(sscanf(str,"%04d%02d%02d/%02d%02d %d %lf %lf %d", 
		     &yr,&mnt,&day,&hour,&min, 
		     &site, &meanval, &sd, &n)==9)
	{
	  HydDateTime dt(yr,mnt,day,hour,min),start=dt.StartOfYear();
	  long int secs=dt-start;
	  double tm=(double)yr + double(secs)/
	    double(dt.daysInYear())/1440.0/60.0;

	  if(use_site[site] && (start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    // increase the number of lines to be read
	    lines++;
	}
      else if(sscanf(str,"%04d%02d%02d/%02d%02d %d %lf %lf", 
		     &yr,&mnt,&day,&hour,&min, 
		     &site, &meanval, &sd)==8)
	{
	  HydDateTime dt(yr,mnt,day,hour,min),start=dt.StartOfYear();
	  long int secs=dt-start;
	  double tm=(double)yr + double(secs)/
	    double(dt.daysInYear())/1440.0/60.0;

	  if(use_site[site] && (start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    // increase the number of lines to be read
	    lines++;
	}
      else if(sscanf(str,"%04d%02d%02d/%02d%02d %d %lf", 
		     &yr,&mnt,&day,&hour,&min, 
		     &site, &meanval)==7 && count_spaces(str)==2)
	{
	  HydDateTime dt(yr,mnt,day,hour,min),start=dt.StartOfYear();
	  long int secs=dt-start;
	  double tm=(double)yr + double(secs)/
	    double(dt.daysInYear())/1440.0/60.0;

	  if(use_site[site] && (start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    // increase the number of lines to be read
	    lines++;
	}
      else if(sscanf(str,"%04d%02d%02d/%02d%02d %lf", &yr,&mnt,&day,&hour,&min, 
		     &meanval)==6)
	{
	  HydDateTime dt(yr,mnt,day,hour,min),start=dt.StartOfYear();
	  long int secs=dt-start;
	  double tm=(double)yr + double(secs)/
	    double(dt.daysInYear())/1440.0/60.0;

	  if((start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    // increase the number of lines to be read
	    lines++;
	}
      
      in.getline(str,9999);
    }
  in.close(); // close the file
  
  // if no lines were found, tell the user and exit:
  if(lines==0)
    {
#ifdef MAIN
      cerr << "No lines read. Illegal format?" << endl;
      exit(0);
#else
      Rcout << "No lines read. Illegal format?" << std::endl;
#endif // MAIN
    }
  
  // Make the measurement array:
  *len=lines;
  series_measurements *ret=new series_measurements[lines];
  
  in.open(infile, ios::in);
  i=0; // counter
  if(in.fail())
    {
#ifdef MAIN
      printf("Couldn't open file \"%s\"", infile);
      exit(0);
#else
      Rcout << "Couldn't open file \"" << infile << "\"!" << std::endl;
#endif // MAIN
    }
  
  in.getline(str,9999);
  while(!in.eof())
    {
      // attempt to read the contents:
      if(sscanf(str,"%04d-%02d-%02d %02d:%02d:%02d;%lf", 
		     &yr,&mnt,&day,&hour,&min,&sec, 
		     &meanval)==7)
	{
	  HydDateTime dt(yr,mnt,day,hour,min,sec),start=dt.StartOfYear(),ref(1970);
	  int secs=dt-start;
	  double tm=(double)yr + double(secs)/
	    double(dt.daysInYear())/1440.0/60.0;
	  
	  if((start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    {
	      // put the contents into the current measurement element:
	      ret[i].site=0;
	      ret[i].tm=double(dt-ref);
	      ret[i].dt=dt;
	      ret[i].meanval=meanval;
	      ret[i].sd=MISSING_VALUE;
	      ret[i].n=(int) MISSING_VALUE;

	      i++; // update the array index
	    }
	    
	}
      else if(!nosites && sscanf(str,"%d %lf %lf %lf %d", &site, 
		&tm, &meanval,&sd,&n)==5)
	{
	  // if this went well and the measurement is not to be removed...
	  if(use_site[site] && (start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    {
	      // put the contents into the current measurement element:
	      ret[i].site=site;
	      ret[i].tm=tm;
	      ret[i].dt=NoHydDateTime;
	      ret[i].meanval=meanval;
	      ret[i].sd=sd;
	      ret[i].n=n;
	      i++; // update the array index
	    }
	}
      else if(!nosites && sscanf(str,"%d %lf %lf %lf", &site, 
		     &tm, &meanval,&sd)==4)
	{
	  if(use_site[site] && (start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    {
	      // put the contents into the current measurement element:
	      ret[i].site=site;
	      ret[i].tm=tm;
	      ret[i].dt=NoHydDateTime;
	      ret[i].meanval=meanval;
	      ret[i].sd=sd;
	      ret[i].n=(int) MISSING_VALUE;
	      i++; // update the array index
	    }
	}
      else if(!nosites && sscanf(str,"%d %lf %lf", &site,&tm, &meanval)==3 &&
	      getnextint(str,&site)[0]==' ')
	{
	  if(use_site[site] && (start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    {
	      // put the contents into the current measurement element:
	      ret[i].site=site;
	      ret[i].tm=tm;
	      ret[i].dt=NoHydDateTime;
	      ret[i].meanval=meanval;
	      ret[i].sd=MISSING_VALUE;
	      ret[i].n=(int) MISSING_VALUE;
	      i++; // update the array index
	    }
	}
      else if(sscanf(str,"%lf %lf %lf %d", &tm, &meanval, &sd, &n)==4)
	{
	  if((start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    {
	      // put the contents into the current measurement element:
	      ret[i].site=0;
	      ret[i].tm=tm;
	      ret[i].dt=NoHydDateTime;
	      ret[i].meanval=meanval;
	      ret[i].sd=sd;
	      ret[i].n=n;
	      i++; // update the array index
	    }
	}
      else if(sscanf(str,"%lf %lf %lf", &tm, &meanval, &sd)==3)
	{
	  if((start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    {
	      // put the contents into the current measurement element:
	      ret[i].site=0;
	      ret[i].tm=tm;
	      ret[i].dt=NoHydDateTime;
	      ret[i].meanval=meanval;
	      ret[i].sd=sd;
	      ret[i].n=(int) MISSING_VALUE;
	      i++; // update the array index
	    }
	}
      else if(sscanf(str,"%lf %lf", &tm, &meanval)==2)
	{
	  if((start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    {
	      // put the contents into the current measurement element:
	      ret[i].site=0;
	      ret[i].tm=tm;
	      ret[i].dt=NoHydDateTime;
	      ret[i].meanval=meanval;
	      ret[i].sd=MISSING_VALUE;
	      ret[i].n=(int) MISSING_VALUE;
	      i++; // update the array index
	    }
	}
      else if(sscanf(str,"%04d%02d%02d/%02d%02d %d %lf %lf %d", 
		     &yr,&mnt,&day,&hour,&min, 
		     &site, &meanval, &sd, &n)==9)
	{
	  HydDateTime dt(yr,mnt,day,hour,min),start=dt.StartOfYear();
	  int secs=dt-start;
	  double tm=
	    (double)yr + double(secs)/double(dt.daysInYear())/1440.0/60.0;
	  
	  if(use_site[site] && (start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    {
	      // put the contents into the current measurement element:
	      ret[i].site=site;
	      ret[i].tm=secs;
	      ret[i].dt=dt;
	      ret[i].meanval=meanval;
	      ret[i].sd=sd;
	      ret[i].n=n;
	      i++; // update the array index
	    }
	}
      else if(sscanf(str,"%04d%02d%02d/%02d%02d %d %lf %lf", 
		     &yr,&mnt,&day,&hour,&min, 
		     &site, &meanval, &sd)==8)
	{
	  HydDateTime dt(yr,mnt,day,hour,min),start=dt.StartOfYear(),ref(1970);
	  int secs=dt-start;
	  double tm=(double)yr + double(secs)/
	    double(dt.daysInYear())/1440.0/60.0;

	  if(use_site[site] && (start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    {
	      // put the contents into the current measurement element:
	      ret[i].site=site;
	      ret[i].tm=double(dt-ref);
	      ret[i].dt=dt;
	      ret[i].meanval=meanval;
	      ret[i].sd=sd;
	      ret[i].n=(int) MISSING_VALUE;
	      i++; // update the array index
	    }
	}
      else if(sscanf(str,"%04d%02d%02d/%02d%02d %d %lf", 
		     &yr,&mnt,&day,&hour,&min, 
		     &site, &meanval)==7 && count_spaces(str)==2)
	{
	  HydDateTime dt(yr,mnt,day,hour,min),start=dt.StartOfYear(),ref(1970);
	  int secs=dt-start;
	  double tm=(double)yr + double(secs)/
	    double(dt.daysInYear())/1440.0/60.0;

	  if(use_site[site] && (start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    {
	      // put the contents into the current measurement element:
	      ret[i].site=site;
	      ret[i].tm=double(dt-ref);
	      ret[i].dt=dt;
	      ret[i].meanval=meanval;
	      ret[i].sd=MISSING_VALUE;
	      ret[i].n=(int) MISSING_VALUE;
	      i++; // update the array index
	    }
	}
      else if(sscanf(str,"%04d%02d%02d/%02d%02d %lf", &yr,&mnt,&day,&hour,&min, 
		     &meanval)==6)
	{
	  HydDateTime dt(yr,mnt,day,hour,min),start=dt.StartOfYear(),ref(1970);
	  int secs=dt-start;
	  double tm=(double)yr + double(secs)/
	    double(dt.daysInYear())/1440.0/60.0;

	  if((start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	    {
	      // put the contents into the current measurement element:
	      ret[i].site=0;
	      ret[i].tm=double(dt-ref);
	      ret[i].dt=dt;
	      ret[i].meanval=meanval;
	      ret[i].sd=MISSING_VALUE;
	      ret[i].n=(int) MISSING_VALUE;

	      i++; // update the array index
	    }
	    
	}

      in.getline(str,9999);
    }
  in.close();

  if(first)
    {
      // count the number of sites:
      numsites=1;
      for(i=0;i<(*len);i++)
	if((int)numsites<(ret[i].site+1))
	  numsites=ret[i].site+1;
      first=0;

      if(!silent)
#ifdef MAIN
	cout << "Numsites=" << numsites << endl; 
#else
	Rcout << "Numsites=" << numsites << std::endl; 
#endif // MAIN

      if(some_pairwise_correlations && numsites>2)
	{
	  // Estimate the probability of a correlation drawn from the prior
	  // resulting in a positively definite covariance matrix.
	  // Used for correcting this prior to this restriction.
	  /* long int t1=clock();
	  p_pos_sigma2=p_positive_correlations(numsites,1000000);
	  cout << "p_pos_sigma2=" << p_pos_sigma2 << endl;
	  long int t2=clock();
	  cout << "Computing time=" << double(t2-t1)/
	  double(CLOCKS_PER_SEC) << endl; */
#ifdef MAIN
	  long int t3=clock();
	  p_pos_site_sigma2=p_positive_sitecorr(numsites,10000);
	  cout << "p_pos_site_sigma2=" << p_pos_site_sigma2 << endl;
	  long int t4=clock();
	  cout << "Computing time=" << double(t4-t3)/
	    double(CLOCKS_PER_SEC) << endl;
	  //exit(0);
#else 
	  p_pos_site_sigma2=p_positive_sitecorr(numsites,10000);
#endif // MAIN
	}
      else
	p_pos_site_sigma2=1.0;
    }
  
  // return the array:
  return ret;
}



// get_correlations: read the series correlation file and store it in
// a timeseries of arrays, corresponding to corr_1_2, corr_1_3, ... ,
// coor_1_S, corr_2_3, ..., corr_2_S, ...,  corr_S-1_S. In total S*(S-1)/2
// correlations, where S is the number of series. Allows for removing some 
// measurements, by site or by age.
double **get_correlations(char *infile, int *len, int **sites, 
			  double **t, bool datetime_format, 
			  bool *use_site,
			  double start_time, 
			  double end_time)
{
  int S=num_series, array_size=S*(S-1)/2;
  // buffer variables for the contents of each line:
  double tm,corr; 
  int site;
  int yr,mnt,day,hour,min; // only for recent time series
  // (time will be counted in minutes rather than m,illions of years)
  
  int i,lines=0; // index variables

  //Rcout << "get_correlations" << std::endl;
  
  ifstream in;
  char str[10000];

  in.open(infile, ios::in);
  if(in.fail())
    {
#ifdef MAIN
      printf("Couldn't open file \"%s\"", infile);
      exit(0);
#else
      Rcout << "Couldn't open file \"" << infile << "\"!" << std::endl;
#endif // MAIN
    }
  
  in.getline(str,9999);
  while(!in.eof())
    {
      char *ptr=str, sbuffer[100];

      if(numsites>1)
	ptr=getnextint(ptr, &site);

      if(ptr)
	{
	  if(datetime_format)
	    {
	      ptr=getnextstr(ptr,sbuffer);
	      if(sscanf(sbuffer,"%04d%02d%02d/%02d%02d", 
			&yr,&mnt,&day,&hour,&min)!=5)
		ptr=NULL; 
	      else
		{
		  HydDateTime dt(yr,mnt,day,hour,min),start=dt.StartOfYear();
		  if(!dt.legal())
		    ptr=NULL;
		  else
		    {
		      long int secs=dt-start;
		      tm=(double)yr + double(secs)/double(dt.daysInYear())/1440.0/60.0;
		    }
		}
	    }
	  else
	    ptr=getnextdouble(ptr,&tm);

	  if(ptr)
	    for(i=0;i<array_size && ptr!=NULL;i++)
	      ptr=getnextdouble(ptr,&corr);
	      
	}

      if(ptr!=NULL && (numsites==1 || use_site[site]) && 
	 (start_time==MISSING_VALUE || tm>start_time) &&
	 (end_time==MISSING_VALUE || tm<end_time))
	lines++;

      in.getline(str,9999);
    }
  in.close(); // close the file
  
  // if no lines were found, tell the user and exit:
  if(lines==0)
    {
#ifdef MAIN
      printf("No lines read in the correlation file. Illegal format?\n");
      exit(0);
#else
      Rcout << "No lines read in the correlation file. Illegal format?"
	    << std::endl;
#endif // MAIN
    }
  
  // Make the correlation timeseries array of arrays:
  double **ret=Make_matrix(lines,array_size);
  double *T=new double[lines];
  int *s=new int[lines];
  
  // Read the contents into the measurement array:
  in.open(infile, ios::in);
  int j=0; // counter
  if(in.fail())
    {
#ifdef MAIN
      printf("Couldn't open file \"%s\"", infile);
      exit(0);
#else
      Rcout << "Couldn't open file \"" << infile << "\"!" << std::endl;
#endif // MAIN
    }
  
  in.getline(str,9999);
  while(!in.eof())
    {
      char *ptr=str, sbuffer[100];
      double *corrbuffer=new double[array_size];

      if(numsites>1)
	ptr=getnextint(ptr, &site);
      
      if(ptr)
	{
	  if(datetime_format)
	    {
	      ptr=getnextstr(ptr,sbuffer);
	      if(sscanf(sbuffer,"%04d%02d%02d/%02d%02d", 
			&yr,&mnt,&day,&hour,&min)!=5)
		ptr=NULL; 
	      else
		{
		  HydDateTime dt(yr,mnt,day,hour,min),start=dt.StartOfYear();
		  if(!dt.legal())
		    ptr=NULL;
		  else
		    {
		      long int secs=dt-start;
		      tm=(double)yr + double(secs)/double(dt.daysInYear())/1440.0/60.0;
		    }
		}
	    }
	  else
	    ptr=getnextdouble(ptr,&tm);
	  
	  if(ptr)
	    for(i=0;i<array_size && ptr!=NULL;i++)
	      ptr=getnextdouble(ptr,&(corrbuffer[i]));
	      
	}

      if(ptr!=NULL && (numsites==1 || use_site[site]) && 
	 (start_time==MISSING_VALUE || tm>start_time) &&
	     (end_time==MISSING_VALUE || tm<end_time))
	{
	  T[j]=tm;
	  if(numsites>1)
	    s[j]=site;
	  else
	    s[j]=0;

	  for(i=0;i<array_size;i++)
	    ret[j][i]=corrbuffer[i];
	  
	  // Test if the correlation matrix is positively definite
	  double **corr_matrix=Make_matrix(S,S);
	  
	  for(i=0;i<S;i++)
	    corr_matrix[i][i]=1.0;

	  i=0;
	  for(int k=0;k<(S-1);k++)
	    for(int l=i+1;l<S;l++)
	      {
		corr_matrix[k][l]=corr_matrix[l][k]=corrbuffer[i];
		i++;
	      }
	  
	  double *corr_lambda=double_eigenvalues(corr_matrix, S);
			
	  for(i=0;i<S;i++)
	    if(corr_lambda[i]<0.0 || !(corr_lambda[i]<1e+200))
	      {
#ifdef MAIN
		cerr << "Non-positively definite correlation "
		  "matrix at line \"" << str << "\"" << endl;
		exit(0);
#else
		Rcout << "Non-positively definite correlation "
		  "matrix at line \"" << str << "\"" << std::endl;
#endif // MAIN
	      }
	  
	  doubledelete(corr_matrix,S);
	  delete [] corr_lambda;
	  j++;
	}

      delete [] corrbuffer;
      in.getline(str,9999);
    }
  in.close();

  *t=T;
  *sites=s;
  *len=lines;

  // return the array:
  return ret;
}

// ********************************************************
// transform_paramete: Transforming from original 
// parametrization to transformed parameters.
// ********************************************************
double transform_parameter(double val, transform_type type)
{
  double minus_mincorr=numsites>1 ? 1.0/double(numsites-1) : 1.0;

  if(val==MISSING_VALUE)
    return MISSING_VALUE;
  
  switch(type)
    {
    case T_LIN:
      return val;
      break;
    case T_LOG:
      return log(val);
      break;
    case T_LOG_INV:
      return log(val);
      break;
    case T_INV:
      return 1.0/val;
      break;
    case T_LOGIST_GLOBAL:
      return log((minus_mincorr+val)/(1.0-val));
      break;
    case T_LOGIST_PAIR:
      return log((1.0+val)/(1.0-val));
      break;
    case T_BINARY:
      if(val!=0.0)
	return 1.0;
      else
	return 0.0;
      break;
    default:
      return MISSING_VALUE;
      break;
    }
  
  return 0.0;
}


// ********************************************************
// invtransform_paramete: Transforming from transformed
// parametrization to original parameters.
// ********************************************************
double invtransform_parameter(double val, transform_type type)
{
  double minus_mincorr=numsites>1 ? 1.0/double(numsites-1) : 1.0;

  if(val==MISSING_VALUE)
    return MISSING_VALUE;
  
  switch(type)
    {
    case T_LIN:
      return val;
      break;
    case T_LOG:
      return exp(val);
      break;
    case T_LOG_INV:
      return exp(val);
      break;
    case T_INV:
      return 1.0/val;
      break;
    case T_LOGIST_GLOBAL:
      return (exp(val)-minus_mincorr)/(1.0+exp(val));
      break;
    case T_LOGIST_PAIR:
      return (exp(val)-1.0)/(1.0+exp(val));
      break;
    case T_BINARY:
      if(val!=0.0)
	return 1.0;
      else
	return 0.0;
      break;
    default:
      return MISSING_VALUE;
      break;
    }

  return 0.0;
}



void series::init(void)
{
  serie_num=0;
  no_layers=false;
  num_per=0;

  for(int i=0;i<100;i++)
    regional_pull[i]=no_sigma[i]=time_integral[i]=regional_sigma[i]=
      sigma_correlated[i]=sigma_pairwise_correlated[i]=sigma_1dim[i]=
      indicator_corr[i]=indicator_corr2[i]=
      indicator_sigma[i]=indicator_pull[i]=0;
  
  meas=comp_meas=NULL;
  comp_len=meas_len=0;
  
  linear_time_dep=false;
  indicator_mu=indicator_lin_t=0;
  
  init_treatment=0;
  init_time=MISSING_VALUE;
  init_dt=NoHydDateTime;
  regional_init=1;
  layered_init=1;
  init_value=MISSING_VALUE;
  
  numlayers=1;
  no_pull_lower=0;

  pr=NULL;
  
  strcpy(name,"X");

  mu=corr=lin_t=NULL;
  init_mu=pull=sigma=NULL;
  paircorr=NULL;
  
  regional_mu=regional_lin_t=0;
  indicator_mu=indicator_lin_t=0;
  allow_positive_pulls=false;
  beta=0.0;
  obs_sd=1.0;
}

void series::cleanup(void)
{
  if(meas)
    delete [] meas;
  if(comp_meas)
    delete [] comp_meas;
  if(pr)
    delete pr;
  if(mu)
    delete [] mu;
  if(corr)
    delete [] corr;
  if(lin_t)
    delete [] lin_t;
  
  doubledelete(init_mu,this->numlayers);
  doubledelete(pull,this->numlayers);
  doubledelete(sigma,this->numlayers);
  tripledelete(paircorr,this->numlayers,numsites);
  
  init();
}

series::series()
{
  init();
}

series::~series()
{
  cleanup();
}

void series::set_series(char *serie_name, int serienum,
			double *X_time, double *X_value,
			unsigned int n, unsigned int num_layers,
			prior *useprior, bool is_datetime,
			double *sd, int *num_meas_per_value, int *sites)
{
  unsigned int i;
  
  init();

  strcpy(name,serie_name);
  numlayers=num_layers;
  if(num_layers==0)
    {
      numlayers=1;
      no_layers=true;
    }
  
  serie_num=serienum;
  
  meas_len=n;
  meas=new series_measurements[n];
  for(i=0;i<n;i++)
    {
      meas[i].serie_num=serienum;
      if(sites)
	meas[i].site=sites[i];
      else
	meas[i].site=0;
      if(num_meas_per_value)
	meas[i].n=num_meas_per_value[i];
      else
	meas[i].n=1;
      meas[i].tm=X_time[i];
      if(!is_datetime)
	meas[i].dt=NoHydDateTime;
      else
	{
	  HydDateTime dt(1970,1,1,0,0);
	  dt+=(int)X_time[i];
	  meas[i].dt=dt;
	}
      meas[i].meanval=X_value[i];
      if(sd)
	meas[i].sd=sd[i];
      else
	meas[i].sd=MISSING_VALUE;
      meas[i].index=0;
    }
  
  mean_time=0.0;
  mean_val=0.0;
  int k=0;
  for(i=0;i<meas_len;i++)
    {
      mean_time+=meas[i].tm;
      if(meas[i].meanval!=MISSING_VALUE)
	{
	  mean_val+=meas[i].meanval;
	  k++;
	}
    }
  mean_time/=double(meas_len);
  mean_val/=double(k);
  
  if(useprior)
    pr=new prior(useprior);
  else
    pr=new prior(-10.0,10.0, 0.001,1000.0, 0.01,10.0, -1.0,1.0, -1.0,1.0, 
		 -10.0,10.0, 0.01,1.0, 1.0);
  
  if(!silent)
    pr->show();
}

void series::set_series(double **X, unsigned int n, unsigned int num_layers)
{
  unsigned int i;
  
  init();

  numlayers=num_layers;
  if(num_layers==0)
    {
      numlayers=1;
      no_layers=true;
    }
  
  meas_len=n;
  meas=new series_measurements[n];
  for(i=0;i<n;i++)
    {
      meas[i].serie_num=0;
      meas[i].site=0;
      meas[i].n=1;
      meas[i].tm=X[i][0];
      meas[i].dt=NoHydDateTime;
      meas[i].meanval=X[i][1];
      meas[i].sd=MISSING_VALUE;
      meas[i].index=0;
    }
  
  mean_time=0.0;
  mean_val=0.0;
  int k=0;
  for(i=0;i<meas_len;i++)
    {
      mean_time+=meas[i].tm;
      if(meas[i].meanval!=MISSING_VALUE)
	{
	  mean_val+=meas[i].meanval;
	  k++;
	}
    }
  mean_time/=double(meas_len);
  mean_val/=double(k);
  
  pr=new prior(-10.0,10.0, 0.001,1000.0, 0.01,10.0, -1.0,1.0, -1.0,1.0, 
	       -10.0,10.0, 0.01,1.0, 1.0);
}


#ifdef MAIN
void series::read_series(char **&argv, int &argc, bool is_main_series,  
			 int serie_number,
			 bool *use_site,double start_time,double end_time)
{
  unsigned int i;
  
  init();
  
  serie_num=serie_number;
  
  // traverse the user options:
  while(argc>1 && argv[1][0]=='-')
    {
      // condition on the second letter:
      switch(argv[1][1])
	{
	case 'P':
	  Tper[num_per]=atof(argv[2]);
	  num_per++;
	  argc--;
	  argv++;
	  break;
	case 't':
	  linear_time_dep=true;
	  break;
	case 'T':
	  {
	    if(argc<3)
	      {
		cerr << "Time integral layer not given!" << endl;
		usage();
	      }
	    
	    unsigned int layer=(unsigned int) (atoi(argv[2])-1);
	    if(layer>=(numlayers-1))
	      {
		cerr << "Can't put time integration on the "
		  "lowest layer (or lower)" << endl;
		exit(0);
	      }
	    time_integral[layer]=1;
	    id_strategy=ID_NONE;
	    no_sigma[layer]=1;
	    argc--;
	    argv++;
	    break;
	  }
	case 'l':
	  if(argc<3)
	    {
	      cerr << "Number of layers not given!" << endl;
	      usage();
	    }
	  numlayers=atoi(argv[2]);
	  if(numlayers==0)
	    {
	      no_layers=true;
	      numlayers=1;
	    }
	  else if(numlayers<0)
	    {
	      cerr << "Not possible to have a negative "
		"number of layers!" << endl;
	      exit(0);
	    }
	  else if(numlayers>=100)
	    {
	      cerr << "A hundred or more layers not possible" << std::endl;
	      exit(0);
	    }
	  argc--;
	  argv++;
	  break;
	case 'c':
	  {
	    comp_meas=get_measurements(argv[2],&comp_len, 
				       use_site, start_time, end_time);
	    argc--;
	    argv++;
	    break;
	  }
	case 'i':
	  {
	    switch(argv[1][2])
	      {
	      case '0':
		init_treatment=1;
		break;
	      case 't':
		init_treatment=1;
		init_time=atof(argv[2]);
		argc--;
		argv++;
		break;
	      case 'T':
		{
		  HydDateTime dt(argv[2]);
		  if(!dt.legal())
		    {
		      cerr << "Illegal date/time format!" << endl;
		      exit(0);
		    }
		  
		  init_treatment=1;
		  init_dt=dt;
		  argc--;
		  argv++;
		  break;
		}
	      case 's':
		regional_init=0;
		break;
	      case 'l':
		layered_init=0;
		break;
	      default:
		cerr << "Unknown initial state option!" << endl;
		usage();
		break;
	      }
	    break;
	  }
	case 'I':
	  {
	    switch(argv[1][2])
	      {
	      case 't':
		init_treatment=1;
		init_time=atof(argv[2]);
		init_value=atof(argv[3]);
		regional_init=layered_init=0;
		argc-=2;
		argv+=2;
		break;
	      case 'T':
		{
		  HydDateTime dt(argv[2]);
		  if(!dt.legal())
		    {
		      cerr << "Illegal date/time format!" << endl;
		      exit(0);
		    }

		  init_treatment=1;
		  init_dt=dt;
		  init_value=atof(argv[3]);
		  regional_init=layered_init=0;
		  argc-=2;
		  argv+=2;
		  break;
		}
		break;
	      default:
		cerr << "Unknown initial state option!" << endl;
		usage();
		break;
	      }
	    break;
	  }
	case 'U':
	  if(!init_treatment)
	    {
	      cerr << "Positive pulls allowed without having any "
		"initial value treatment!" << endl;
	      exit(0);
	    }
	  allow_positive_pulls=true;
	  id_strategy=ID_NONE;
	  break;
	case 'S':
	  {
	    int layer;
	    switch(argv[1][2])
	      {
	      case 'S':
		layer=atoi(argv[2])-1;
		indicator_corr[layer]=1;
		argc--;
		argv++;
		break;
	      case 'C':
		layer=atoi(argv[2])-1;
		indicator_corr2[layer]=1;
		argc--;
		argv++;
		break;
	      case 's':
		layer=atoi(argv[2])-1;
		indicator_sigma[layer]=1;
		argc--;
		argv++;
		break;
	      case 'p':
		layer=atoi(argv[2])-1;
		indicator_pull[layer]=1;
		argc--;
		argv++;
		break;
	      case 'u':
		indicator_mu=1;
		break;
	      case 't':
		indicator_lin_t=1;
		linear_time_dep=true;
		break;
	      default:
		cerr << "Unknown indicator option!" << endl;
		usage();
		break;
	      }
	    
	    use_indicator=1;
	    break;
	  }
	case 'r':
	  {
	    switch(argv[1][2])
	      {
	      case 'u':
		regional_mu=1;
		break;
	      case 't':
		regional_lin_t=1;
		linear_time_dep=true;
		break;
	      case 'p':
		{
		  int layer=atoi(argv[2])-1;
		  regional_pull[layer]=1;
		  argc--;
		  argv++;
		  break;
		}
	      case 's':
		{
		  int layer=atoi(argv[2])-1;
		  regional_sigma[layer]=1;
		  argc--;
		  argv++;
		  break;
		}
	      case 'C':
		{
		  int layer=atoi(argv[2])-1;
		  sigma_pairwise_correlated[layer]=1;
		  some_pairwise_correlations=true;
		  argc--;
		  argv++;
		  break;
		}
	      default:
		cerr << "Unknown regional option!" << endl;
		usage();
		break;
	      }

	    break;
	  }
	case 'C':
	  {
	    int layer=atoi(argv[2])-1;
	    sigma_correlated[layer]=1;
	    argc--;
	    argv++;
	    break;
	  }
	case 'n':
	  {
	    switch(argv[1][2])
	      {
	      case 'p':
		no_pull_lower=1;
		break;
	      case 's':
		{
		  int layer=atoi(argv[2])-1;
		  no_sigma[layer]=1;
		  argc--;
		  argv++;
		  break;
		}
	      default:
		cerr << "Unknown no noise option!" << endl;
		usage();
		break;
	      }
	    break;
	  }
	case '1':
	  {
	    switch(argv[1][2])
	      {
	      case 's':
		{
		  int layer=atoi(argv[2])-1;
		  sigma_1dim[layer]=1;
		  argc--;
		  argv++;
		  break;
		}
	      default:
		cerr << "Unknown 1D option!" << endl;
		usage();
		break;
	      }
	    
	    break;
	  }
	default:
	  cerr << "Unknown series option!" << endl;
	  usage();
	  break;
	}

      argc--;
      argv++;
    }
  
  strcpy(name,argv[1]);
  meas=get_measurements(argv[2],&meas_len, use_site, start_time, end_time);

  if(init_treatment && (init_time!=MISSING_VALUE || init_dt!=NoHydDateTime))
    {
      if(init_dt!=NoHydDateTime)
	{
	  HydDateTime start=init_dt.StartOfYear(),ref(1970);
	  // int secs=init_dt-start, yr=init_dt.getYear();;
	  //double tm=(double)yr + double(secs)/double(init_dt.daysInYear())/1440.0/60.0;

	  init_time=double(start-ref);
	}

      if(init_time >= meas[0].tm)
	{
	  cerr << "Initial time comes after first measurement time!" << endl;
	  exit(0);
	}

      series_measurements ms0;
      ms0.site=0;
      ms0.n=ms0.serie_num=ms0.index=0;
      ms0.dt=init_dt;
      ms0.tm=init_time;
      ms0.meanval=ms0.sd=MISSING_VALUE;

      int measlen2=meas_len+1;
      series_measurements *meas2=new series_measurements[measlen2];
      meas2[0]=ms0;
      for(i=0;i<meas_len;i++)
	meas2[i+1]=meas[i];
      
      delete [] meas;
      meas=meas2;
      meas_len=measlen2;
    }

  for(i=0;i<meas_len;i++)
    meas[i].serie_num=serie_number;

  mean_time=0.0;
  mean_val=0.0;
  int k=0;
  for(i=0;i<meas_len;i++)
    {
      mean_time+=meas[i].tm;
      if(meas[i].meanval!=MISSING_VALUE)
	{
	  mean_val+=meas[i].meanval;
	  k++;
	}
    }
  mean_time/=double(meas_len);
  mean_val/=double(k);

  pr=new prior(argv[3],mean_val);
  if(pr->is_log==1)
    {
      for(i=0;i<meas_len;i++)
	{
	  if(meas[i].sd!=MISSING_VALUE)
	    {
	      cerr << "Standard deviations given for a series can  not"
		" apply for the log-transformed series!" << endl;
	      exit(0);
	    }

	  if(meas[i].meanval!=MISSING_VALUE)
	    meas[i].meanval=log(meas[i].meanval);
	}
    }
  
  if(init_treatment && init_value!=MISSING_VALUE)
    {
      pr->init_m=init_value;
      pr->init_s=1e-100;
      pr->init_1=pr->init_m-1.96*pr->init_s;
      pr->init_2=pr->init_m+1.96*pr->init_s;
    }
  
  argc-=2;
  argv+=2; // this is also done inside the main option loop,
  // which is why it's not -3.
}
#endif // MAIN



double ml_indep(double *residuals, int len,double *mu0, double *sd0)
{
  double mu=find_statistics(residuals,len,MEAN);
  double sd=sqrt(double(len-1)/double(len))*
    find_statistics(residuals,len,STANDARD_DEVIATION);
  double n=double(len);
  double ret=-0.5*n*log(2.0*M_PI)-n*log(sd);

  for(int i=0;i<len;i++)
    ret-= 0.5*(residuals[i]-mu)*(residuals[i]-mu)/sd/sd;

  *mu0=mu;
  *sd0=sd;

  return ret;
}

double *ar1_res=NULL;
int ar1_len=0;
double minusloglik_ar1(double *pars)
{
  double mu=pars[0];
  double sd=exp(pars[1]);
  double a=-1.0+2.0*exp(pars[2])/(1.0+exp(pars[2]));
  double n=double(ar1_len);
  double ret=-0.5*n*log(2.0*M_PI)-n*log(sd)-0.5*(n-1)*log(1.0-a*a);
  
  ret -= 0.5*(ar1_res[0]-mu)*(ar1_res[0]-mu)/sd/sd;
  
  for(int i=1;i<ar1_len;i++)
    ret -= 0.5*(ar1_res[i]-mu-a*(ar1_res[i-1]-mu))*
      (ar1_res[i]-mu-a*(ar1_res[i-1]-mu))/
      (sd*sd*(1.0-a*a));

  return -ret;
}

double ml_ar1(double *residuals, int len,double *mu1, double *sd1, double *a)
{
  double mu=MISSING_VALUE,sd=MISSING_VALUE,aa=MISSING_VALUE;
  double *pars=new double[3], best_loglik=MISSING_VALUE;
  int i,j;
  
  if(!ar1_res)
    ar1_res=new double[len];
  ar1_len=len;
  for(i=0;i<len;i++)
    ar1_res[i]=residuals[i];
  
  /*
  mu=find_statistics(residuals,len,MEAN);
  sd=sqrt(double(len-1)/double(len))*
    find_statistics(residuals,len,STANDARD_DEVIATION);
  pars[0]=mu;
  pars[1]=log(sd);
  pars[2]=0.0;
  cout << -minusloglik_ar1(pars) << endl;
  */
  
  for(i=0;i<10;i++)
    {
      for(j=0;j<3;j++)
	pars[j]=get_random_gauss();
      
      // do the optimization:
      int maxiter=1000;
      double *pars2=quasi_newton(minusloglik_ar1, 3, pars, 0.001,maxiter,true,true);
      
      if(maxiter<1000)
	{
	  double log_lik=-minusloglik_ar1(pars2);
	  
	  if(best_loglik==MISSING_VALUE || log_lik>best_loglik)
	    {
	      mu=pars2[0];
	      sd=exp(pars2[1]);
	      aa=-1.0+2.0*exp(pars2[2])/(1.0+exp(pars2[2]));
	      best_loglik=log_lik;
	    }
	}
      
      delete [] pars2;
    }
  
  *a=aa;
  *mu1=mu;
  *sd1=sd;
  
  delete [] pars;
  
  return best_loglik;
}

void test_ar1(double *residuals,  int len)
{
  double mu0,sd0;
  double mu1,sd1,a;
  double ml0=ml_indep(residuals,len,&mu0,&sd0);
  double ml1=ml_ar1(residuals,len,&mu1,&sd1,&a);
  
  double dev=2.0*(ml1-ml0);
  double p_value=1.0-chisq_deg1_cdf(dev);
  //double p_value=gsl_cdf_chisq_Q(dev,1.0);

#ifdef MAIN
  printf("Autocorrelation in the data: %f\n",
	 get_auto_correlation(residuals, len));
  printf("Independence: mu=%f sd=%f      ML=%f\n", mu0,sd0, ml0);
  printf("AR(1):        mu=%f sd=%f a=%f ML=%f\n", mu1,sd1, a,ml1);
  printf("Deviance: %f\n", dev);
  printf("p-value : %g\n\n",p_value);
#else
  Rcout << "Autocorrelation in the data: " <<
    get_auto_correlation(residuals, len) << std::endl;
  Rcout << "Independence: mu=" << mu0 << " sd=" << sd0 << "      ML=" <<
    ml0 << std::endl;
  Rcout << "AR(1)         mu=" << mu1 << " sd=" << sd1 << " a=" << a <<
    " ML=" << ml1 << std::endl;
  Rcout << "Deviance: " << dev << std::endl;
  Rcout << "p-value : " << p_value << std::endl << std::endl;
#endif // MAIN
}

double *ou_res=NULL, *ou_tm=NULL;
int ou_len=0;
double minusloglik_ou(double *pars)
{
  double mu=pars[0];
  double sd=exp(pars[1]);
  double dt1=exp(pars[2]);
  double n=double(ou_len);
  double ret=-0.5*n*log(2.0*M_PI);
  
  ret -= log(sd)+0.5*(ou_res[0]-mu)*(ou_res[0]-mu)/sd/sd;
  
  for(int i=1;i<ou_len;i++)
    ret -= log(sd)+0.5*log(1.0-exp(-2.0*(ou_tm[i]-ou_tm[i-1])/dt1))+
      0.5*(ou_res[i]-mu-exp(-(ou_tm[i]-ou_tm[i-1])/dt1)*(ou_res[i-1]-mu))*
      (ou_res[i]-mu-exp(-(ou_tm[i]-ou_tm[i-1])/dt1)*(ou_res[i-1]-mu))/
      (sd*sd*(1.0-exp(-2.0*(ou_tm[i]-ou_tm[i-1])/dt1)));

  return -ret;
}

double ml_OU(double *residuals, double *tm, int len,double *mu1, double *sd1, double *dt1)
{
  double mu=MISSING_VALUE,sd=MISSING_VALUE,dt=MISSING_VALUE;
  double *pars=new double[3], best_loglik=MISSING_VALUE;
  int i,j;

  if(!ou_res)
    ou_res=new double[len];
  if(!ou_tm)
    ou_tm=new double[len];
  ou_len=len;
  for(i=0;i<len;i++)
    {
      ou_res[i]=residuals[i];
      ou_tm[i]=tm[i];
    }
  
  for(i=0;i<10;i++)
    {
      for(j=0;j<3;j++)
	pars[j]=get_random_gauss();

      // do the optimization:
      int maxiter=1000;
      //double *pars2=gsl_optimization_cover(minusloglik_ou, 3, pars, 0.001,maxiter);
      double *pars2=quasi_newton(minusloglik_ou, 3, pars, 0.001,maxiter);

      if(maxiter<1000)
	{
	  double log_lik=-minusloglik_ou(pars2);
	  
	  if(best_loglik==MISSING_VALUE || log_lik>best_loglik)
	    {
	      mu=pars2[0];
	      sd=exp(pars2[1]);
	      dt=exp(pars2[2]);
	      best_loglik=log_lik;
	    }
	}

      delete [] pars2;
    }

  *dt1=dt;
  *mu1=mu;
  *sd1=sd;

  delete [] pars;

  return best_loglik;
}

void test_OU(double *residuals, double *tm,  int len)
{
  
#ifdef DETAILED_TIMERS
  timers[24][0]=clock();
#endif // DETAILED_TIMERS
  
  double mu0,sd0;
  double mu1,sd1,dt1;
  double ml0=ml_indep(residuals,len,&mu0,&sd0);
  double ml1=ml_OU(residuals,tm,len,&mu1,&sd1,&dt1);
  
  double dev=2.0*(ml1-ml0);
  double p_value=1.0-chisq_deg1_cdf(dev);
  //double p_value=gsl_cdf_chisq_Q(dev,1.0);

#ifdef MAIN
  printf("OU vs independence:\n");
  printf("Independence: mu=%f sd=%f      ML=%f\n", mu0,sd0, ml0);
  printf("AR(1):        mu=%f sd=%f dt1=%f ML=%f\n", mu1,sd1, dt1,ml1);
  printf("Deviance: %f\n", dev);
  printf("p-value : %g\n\n",p_value);
#else
  Rcout << "OU vs independence:" << std::endl;
  Rcout << "Independence: mu=" << mu0 << " sd=" << sd0 << "      ML="
	<< ml0 << std::endl;
  Rcout << "AR(1):        mu=" << mu1 << " sd=" << sd1 << " dt1=" << dt1 <<
    " ML=" << ml1 << std::endl;
  Rcout << "Deviance: " << dev << std::endl;
  Rcout << "p-value : " << p_value << std::endl;
#endif // MAIN
  
#ifdef DETAILED_TIMERS
  timers[24][1]=clock();
  timers[24][2]+=(timers[24][1]-timers[24][0]);
#endif // DETAILED_TIMERS
  
}

void autocorr_analyzer(double *residuals, int len, 
		       int max_corr_len, char *filestart, bool absvalues)
{
  double mu=find_statistics(residuals,len,MEAN);
  double sd=find_statistics(residuals,len,VARIATION);
  double *autocorr=new double[max_corr_len+1];
  double corr_limit=1.96/sqrt(double(len));
  int i,k;

  autocorr[0]=MISSING_VALUE;
  for(k=1;k<=max_corr_len;k++)
    {
      autocorr[k]=0.0;
      for(i=0;i<len-k;i++)
	autocorr[k] += (residuals[i]-mu)*(residuals[i+k]-mu);
      autocorr[k]/=sd*sd*double(len-k);
    }

#ifdef MAIN
  char filename[1000];
  if(!silent)
    {
      FILE *p;

      if(!filestart || !*filestart)
	p=popen("vvgraph -x lag -y autocorr","w");
      else
	{
	  if(!absvalues)
	    snprintf(filename,999,"%s_res_autocorr.txt", filestart);
	  else
	    snprintf(filename,999,"%s_absres_autocorr.txt", filestart);
	  p=fopen(filename,"w");
	}
      fprintf(p,"# Column 1: auto-correlation\n");
      fprintf(p,"# Column 1 - type: filled bars\n");
      fprintf(p,"# Column 2: 95%% confidence bands\n");
      fprintf(p,"########################\n");
      for(k=0;k<=max_corr_len;k++)
	fprintf(p,"%d %f -10000000\n", k, autocorr[k]);
      fprintf(p,"%d -10000000 %f\n", 0, corr_limit);
      fprintf(p,"%d -10000000 %f\n", max_corr_len, corr_limit);
      fprintf(p,"%d -10000000 %f\n", 0, MISSING_VALUE);
      fprintf(p,"%d -10000000 %f\n", 0, -corr_limit);
      fprintf(p,"%d -10000000 %f\n", max_corr_len, -corr_limit);
      pclose(p);
    }
#endif // MAIN
  
  double p_value=2.0*(1.0-standard_normal_cdf(ABSVAL((autocorr[1]*sqrt(double(len))))));

#ifdef MAIN
  printf("Autocorrelation in the data: %f\n",
	 autocorr[1]);
  printf("95%% confidence band: +/-%f\n", corr_limit);
  printf("p-value : %g\n\n",p_value);
#else
  Rcout << "Autocorrelation in the data: " << autocorr[1] << std::endl;
  Rcout << "95% confidence band: +/-" << corr_limit << std::endl;
  Rcout << "p-value : " << p_value << std::endl;
#endif // MAIN
}

void runs_test(double *residuals, int len)
{
  
#ifdef DETAILED_TIMERS
  timers[25][0]=clock();
#endif // DETAILED_TIMERS
  
  int i,R1=1,R2=1, l1=0,l2=0;
  double n=double(len);
  double median=find_statistics(residuals, len, MEDIAN);
  int *run1=new int[len], *run2=new int[len-1];

  for(i=0;i<len;i++)
    run1[i]=residuals[i]>median ? 1 : -1;
  for(i=0;i<len;i++)
    if(run1[i]>0)
      l1++;

  for(i=1;i<len;i++)
    if(run1[i]!=run1[i-1])
      R1++;

  for(i=0;i<(len-1);i++)
    run2[i]=residuals[i+1]>residuals[i] ? 1 : -1;      
  for(i=0;i<(len-1);i++)
    if(run2[i]>0)
      l2++;

  for(i=1;i<(len-1);i++)
    if(run2[i]!=run2[i-1])
      R2++;

  double lambda1=double(l1)/n; // lambda2=double(l2)/(n-1.0);
  double R_stat1=(R1-2.0*n*lambda1*(1.0-lambda1))/(2.0*sqrt(n)*lambda1*(1.0-lambda1));
  //double R_stat2=(R2-2.0*(n-1.0)*lambda2*(1.0-lambda2))/
  // (2.0*sqrt(n-1.0)*lambda2*(1.0-lambda2));
  double p_value1=2.0*(1.0-standard_normal_cdf(ABSVAL((R_stat1))));
  //double p_value1=2.0*gsl_cdf_ugaussian_Q(ABSVAL((R_stat1)));
  //double p_value2=2.0*gsl_cdf_ugaussian_Q(ABSVAL((R_stat2)));
  int lower_limit1=(int) 
    floor(-1.96*2.0*sqrt(n)*lambda1*(1.0-lambda1)+2.0*n*lambda1*(1.0-lambda1));
  int upper_limit1=(int) 
    ceil(+1.96*2.0*sqrt(n)*lambda1*(1.0-lambda1)+2.0*n*lambda1*(1.0-lambda1));
  /*
    int lower_limit2=(int) 
    floor(-1.96*2.0*sqrt(n-1.0)*lambda2*(1.0-lambda2)+2.0*(n-1.0)*lambda2*(1.0-lambda2));
  int upper_limit2=(int) 
    ceil(+1.96*2.0*sqrt(n-1.0)*lambda2*(1.0-lambda2)+2.0*(n-1.0)*lambda2*(1.0-lambda2));
  */

#ifdef MAIN
  printf("Runs test - median:\n");
  printf("Runs: %d of %d  statistics=%f\n",R1,len,R_stat1);
  printf("lambda=%f 1-lambda=%f\n", lambda1, 1.0-lambda1);
  printf("95%% confidence limit: runs:%d-%d , statistics=+/-1.96 \n", 
	 lower_limit1,upper_limit1);
  printf("p-value: %g\n\n", p_value1);
  /*
  printf("Runs test - differance:\n");
  printf("Runs: %d of %d  statistics=%f\n",R2,len-1,R_stat2);
  printf("lambda=%f 1-lambda=%f\n", lambda2, 1.0-lambda2);
  printf("95%% confidence limit: runs:%d-%d , statistics=+/-1.96 \n", 
	 lower_limit2,upper_limit2);
  printf("p-value: %g\n\n", p_value2);
  */
#else
  Rcout << "Runs test - median:" << std::endl;
  Rcout << "Runs: " << R1 << " of " << len << "  statistics=" <<
    R_stat1 << std::endl;
  Rcout << "lambda=" << lambda1 << " 1-lambda=" << 1.0-lambda1 << std::endl;
  Rcout << "95% confidence limit: runs:" << lower_limit1 << "-" <<
    upper_limit1 << " , statistics=+/-1.96 " << std::endl; 
  Rcout << "p-value : " << p_value1 << std::endl;
#endif // MAIN
  
#ifdef DETAILED_TIMERS
  timers[25][1]=clock();
  timers[25][2]+=(timers[25][1]-timers[25][0]);
#endif // DETAILED_TIMERS
  
}

void analyze_residuals(double *residuals, double *tm, HydDateTime *dt, int len,
		       char *filestart)
{
#ifdef DETAILED_TIMERS
  timers[26][0]=clock();
#endif // DETAILED_TIMERS
  
  int i;
  double *abs_res=new double[len];

#ifdef MAIN
  printf("*******************\nResiduals:\n*******************\n\n");
#else
  Rcout << "*******************\nResiduals:\n*******************" << std::endl;
#endif // MAIN
  
  for(i=0;i<len;i++)
    abs_res[i]=ABSVAL((residuals[i]));

  double *resid_ordered=new double[len];
  for(i=0;i<len;i++)
    resid_ordered[i]=residuals[i];
  qsort(resid_ordered, size_t(len),sizeof(double),compare_double);

#ifdef MAIN
  FILE *p,*p2;
  char cmd[1000], filename[1000];
  if(!silent)
    {
      if(!filestart || !*filestart)
	{
	  if(dt && dt[0]!=NoHydDateTime)
	    snprintf(cmd, 999,"timeseriegraph -x t -y res");
	  else
	    snprintf(cmd, 999,"vvgraph -x t -y res");
	  p=popen(cmd,"w"); 
	  p2=popen("histogramme","w");
	}
      else
	{
	  snprintf(filename, 999,"%s_residuals.txt",filestart);
	  p=fopen(filename,"w");
	  snprintf(filename, 999,"%s_res_hist.txt",filestart);
	  p2=fopen(filename,"w");
	}
      
      fprintf(p, "# Column 1: residuals\n");
      fprintf(p, "# Column 1 - type: dot\n");
      fprintf(p, "########################\n");
      
      for(i=0;i<len;i++)
	{
	  if(dt && dt[i]!=NoHydDateTime)
	    fprintf(p,"%s %f\n", dt[i].syCh(5), residuals[i]);
	  else
	    fprintf(p,"%f %f\n", tm[i], residuals[i]);
	  fprintf(p2,"%f\n", residuals[i]);
	}

      if(!filestart || !*filestart)
	{
	  pclose(p);
	  pclose(p2);
	}
      else
	{
	  fclose(p);
	  fclose(p2);
	}
      
      if(!filestart || !*filestart)
	{
	  snprintf(cmd, 999,"vvgraph -x \"normal quantiles\" -y \"residuals\"");
	  p=popen(cmd,"w");
	}
      else
	{
	  snprintf(filename, 999,"%s_res_qqplot.txt", filestart);
	  p=fopen(filename,"w");
	}
      fprintf(p,"# Column 1: residuals\n");
      fprintf(p,"# Column 1 -type : dots\n");
      fprintf(p,"# Column 2: line\n");
      fprintf(p, "########################\n");
      for(i=0;i<len;i++)
	{
	  fprintf(p, "%f %f -10000000\n", 
		  normal_invcdf(double(i+1)/double(len+1),0.0,1.0),
		  //gsl_cdf_ugaussian_Pinv(double(i+1)/double(len+1)),
		  resid_ordered[i]);
	  //printf("%f %f\n", gsl_cdf_ugaussian_Pinv(double(i+1)/double(len+1)),
	  //     resid_ordered[i]);
	}
      double min_val=MINIM(1.0-standard_normal_cdf(double(1)/double(len+1)),
			   resid_ordered[0]);
      double max_val=MAXIM(1.0-standard_normal_cdf(double(len)/double(len+1)),
			   resid_ordered[len-1]);
      //double min_val=MINIM(gsl_cdf_ugaussian_Q(double(1)/double(len+1)),
      //resid_ordered[0]);
      //double max_val=MAXIM(gsl_cdf_ugaussian_Q(double(len)/double(len+1)),
      //resid_ordered[len-1]);
      fprintf(p,"%f -10000000 %f\n", min_val,min_val);
      fprintf(p,"%f -10000000 %f\n", max_val,max_val);
      if(!filestart || !*filestart)
	pclose(p);
      else
	fclose(p);
    }
#endif // MAIN
  
  test_ar1(residuals,len);
  autocorr_analyzer(residuals,len,20,filestart,false);
  test_OU(residuals,tm,len);
  runs_test(residuals,len);

#ifdef MAIN
  printf("*******************\nAbsolute value residuals:\n*******************\n\n");
#else
  Rcout << "*******************" << std::endl << "Absolute value residuals:" <<
    std::endl << "*******************" << std::endl;
#endif // MAIN
  
  test_ar1(abs_res,len);
  autocorr_analyzer(abs_res,len,20,filestart,true);
  test_OU(abs_res,tm,len);
  runs_test(abs_res,len);

#ifdef DETAILED_TIMERS
  timers[26][1]=clock();
  timers[26][2]+=(timers[26][1]-timers[26][0]);
#endif // DETAILED_TIMERS
  
  delete [] abs_res;
  delete [] resid_ordered;
}





void keep_x_and_P(double *t_k,double **x_k_s, double ***P_k_s, unsigned int len)
{
  unsigned int i,j,k;

  if(x_k_s_kept)
    doubledelete(x_k_s_kept,len);
  x_k_s_kept=new double*[len];
  len_x_k_s_kept=len;
  size_x_k_s_kept=num_states;

  if(P_k_s_kept)
    {
      for(k=0;k<len;k++)
	doubledelete(P_k_s_kept[k],num_states);
      delete [] P_k_s_kept;
    }
  P_k_s_kept=new double**[len];
  

  for(k=0;k<len;k++)
    {
      x_k_s_kept[k]=new double[num_states];
      P_k_s_kept[k]=Make_matrix(num_states,num_states);

      for(i=0;i<num_states;i++)
	{
	  x_k_s_kept[k][i]=x_k_s[k][i];
	  for(j=0;j<num_states;j++)
	    P_k_s_kept[k][i][j]=P_k_s[k][i][j];
	}
      t_k_smooth[k]=t_k[k];
      
      if(dt_k_smooth)
	dt_k_smooth[k]=meas_smooth[k].dt;
    }
}

void cleanup_x_and_P(int len)
{
  if(x_k_s_kept)
    doubledelete(x_k_s_kept,len_x_k_s_kept);
  x_k_s_kept=NULL;
  
  if(P_k_s_kept)
    tripledelete(P_k_s_kept,len_x_k_s_kept,size_x_k_s_kept);
  P_k_s_kept=NULL;
  
  len_x_k_s_kept=0;
  size_x_k_s_kept=0;
}


double loglik(double *pars, int dosmooth, int do_realize,
	      int residual_analysis, char *res_filestart, 
	      int debug, char **simulation_files,
	      int return_residuals,
	      double **residuals_time,double ***residuals,
	      double ***prior_expected_values, 
	      int *resid_numcolumns, int *resid_len)
{
  
#ifdef DETAILED_TIMERS
  timers[1][0]=clock();
  timers[2][0]=clock();
#endif // DETAILED_TIMERS
  
  unsigned int s,i,j,k,l,t=0, t_0=0;
  numpar=0;
  bool pairwise_wrong=false;

#ifndef MAIN
  R_CheckUserInterrupt();
#endif // MAIN
  
  /*
  static gsl_rng *rptr=NULL;
  if(rptr==NULL)
    {
    rptr=gsl_rng_alloc(gsl_rng_rand48);
    gsl_rng_set(rptr, rand()); 
    }
  */

#ifndef MAIN
  if(detailed_loglik)
    Rcout << "loglik started" << std::endl;
#endif // MAIN
  
  // **************************************************
  // If the parameter value array is missing, that means
  // that the global variables for storing parameter
  // info are to be filled out. First, the arrays
  // are created. 10000 parameters with parameter names 
  // max 100 byes each should be enough....
  // (Or as Bill Gates would say, 640kb RAM must be
  // enough for all purposes.)
  if(!pars)
    {
      if(par_name)
	doubledelete(par_name,LARGE_ENOUGH_ARRAY);
      par_name=new char*[LARGE_ENOUGH_ARRAY];
      if(par_trans_type)
	delete [] par_trans_type;
      par_trans_type=new transform_type[LARGE_ENOUGH_ARRAY];
      if(par_type)
	delete [] par_type;
      par_type=new param_type[LARGE_ENOUGH_ARRAY];
      if(par_layer)
	delete [] par_layer;
      par_layer=new int[LARGE_ENOUGH_ARRAY];
      if(par_region)
	delete [] par_region;
      par_region=new int[LARGE_ENOUGH_ARRAY];
      if(par_series)
	delete [] par_series;
      par_series=new int[LARGE_ENOUGH_ARRAY];
      for(i=0;i<LARGE_ENOUGH_ARRAY;i++)
	par_name[i]=new char[300];
    
      for(s=0;s<num_series;s++)
	if(ser)
	  {
	    if(!ser[s].mu)
	      ser[s].mu=new double[numsites];
	    if(!ser[s].init_mu)
	      ser[s].init_mu=Make_matrix(ser[s].numlayers,numsites);
	    if(!ser[s].lin_t)
	      ser[s].lin_t=new double[numsites];
	    if(!ser[s].pull)
	      ser[s].pull=Make_matrix(ser[s].numlayers,numsites);
	    if(!ser[s].sigma)
	      ser[s].sigma=Make_matrix(ser[s].numlayers,numsites);
	    if(!ser[s].corr)
	      ser[s].corr=new double[ser[s].numlayers];
	    if(!ser[s].paircorr)
	      {
		ser[s].paircorr=new double**[ser[s].numlayers];
		for(l=0;l<ser[s].numlayers;l++)
		  ser[s].paircorr[l]=Make_matrix(numsites,numsites);
	      }
	  }
      
      if(use_indicator)
	indicator_array=new int[numsites];
    }
  
  // **************************************************
  // If indicator variables are used...
  if(use_indicator)
    {
      // make the indicactor arrays and fill
      // them with approrriate contents:
      if(pars)
	for(i=0;i<numsites;i++) // traverse the sites
	  indicator_array[i]=(int) pars[numpar+i];
      else
	for(i=0;i<numsites;i++) // traverse the sites
	  {
	    par_trans_type[i+numpar]=T_BINARY;
	    par_type[i+numpar]=INDICATOR;
	    par_layer[i+numpar]=ser[0].numlayers+1;
	    par_region[i+numpar]=i;
	    par_series[i+numpar]=0;
	    snprintf(par_name[i+numpar], 99, "indicator%d",i);
	  }
      numpar+=numsites; // update the number of parameters
    }
  

  // **************************************************
  // Read the parameters while handling the model type:

#ifndef MAIN
  if(detailed_loglik)
    Rcout << "start reading parameters" << std::endl;
#endif // MAIN
  
  for(s=0;s<num_series;s++)
    {
      if(ser[s].regional_mu) // regional expectancy?
	{
	  for(i=0;i<numsites;i++) // traverse the sites
	    if(pars) // parameter value array is given?
	      ser[s].mu[i]=pars[i+numpar]; // fetch the parameter value
	    else
	      {
		// fill the global parameter arrays with
		// appropriate contents:
		par_trans_type[i+numpar]=T_LIN;
		par_type[i+numpar]=MU;
		par_layer[i+numpar]=ser[s].numlayers+1;
		par_region[i+numpar]=i;
		par_series[i+numpar]=s;
		snprintf(par_name[i+numpar], 250,"mu_%s_%d",ser[s].name,i);
	      }
	  numpar+=numsites; // update the number of parameters
	}
      else if(ser[s].indicator_mu) // indicator-separated expectancy?
	{
	  // indicator always first parameters, if used
	  
	  if(pars) // parameter value array is given?
	    {
	      // fetch the parameter value
	      for(i=0;i<numsites;i++) // traverse the sites
		if(!indicator_array[i])
		  ser[s].mu[i]=pars[numpar];
		else
		  ser[s].mu[i]=pars[numpar+1];
	    }
	  else
	    {
	      // fill the global parameter arrays with
	      // appropriate contents:
	      for(i=0;i<2;i++)
		{
		  par_trans_type[numpar+i]=T_LIN;
		  par_type[numpar+i]=MU;
		  par_layer[numpar+i]=ser[s].numlayers+1;
		  par_region[numpar+i]=-1; // global
		  par_series[numpar+i]=s;
		  snprintf(par_name[numpar+i], 250, "mu_%s_%d", ser[s].name,i);
		}
	    }
	  
	  numpar+=2; // update the number of parameters
	}
      else // global expectancy
	{
	  if(pars) // parameter value array is given?
	    // fetch the parameter value:
	    for(i=0;i<numsites;i++) // traverse the sites
	      ser[s].mu[i]=pars[numpar];
	  if(!pars)
	    {
	      // fill the global parameter arrays with
	      // appropriate contents:
	      par_trans_type[numpar]=T_LIN;
	      par_type[numpar]=MU;
	      par_layer[numpar]=ser[s].numlayers+1;
	      par_region[numpar]=-1;
	      par_series[numpar]=s;
	      snprintf(par_name[numpar], 250, "mu_%s",ser[s].name);
	    }
      
	  numpar++; // update the number of parameters
	}

      
      if(ser[s].linear_time_dep) // linear time dependency?
	{
	  if(ser[s].regional_lin_t) // regional linear time dependency?
	    {
	      for(i=0;i<numsites;i++) // traverse the sites
		if(pars) // parameter value array is given?
		  ser[s].lin_t[i]=pars[i+numpar]; // fetch the parameter value
		else
		  {
		    // fill the global parameter arrays with
		    // appropriate contents:
		    par_trans_type[i+numpar]=T_LIN;
		    par_type[i+numpar]=LIN_T;
		    par_layer[i+numpar]=ser[s].numlayers+1;
		    par_region[i+numpar]=i;
		    par_series[i+numpar]=s;
		    snprintf(par_name[i+numpar], 250,
			     "lin_t_%s_%d",ser[s].name,i);
		  }
	      numpar+=numsites; // update the number of parameters
	    }
	  else if(ser[s].indicator_lin_t) // indicator-separated linear time dep.?
	    {
	      // indicator always first parameters, if used
	      
	      if(pars) // parameter value array is given?
		{
		  // fetch the parameter value
		  for(i=0;i<numsites;i++) // traverse the sites
		    if(!indicator_array[i])
		      ser[s].lin_t[i]=pars[numpar];
		    else
		      ser[s].lin_t[i]=pars[numpar+1];
		}
	      else
		{
		  // fill the global parameter arrays with
		  // appropriate contents:
		  for(i=0;i<2;i++)
		    {
		      par_trans_type[numpar+i]=T_LIN;
		      par_type[numpar+i]=LIN_T;
		      par_layer[numpar+i]=ser[s].numlayers+1;
		      par_region[numpar+i]=-1;
		      par_series[numpar+i]=s;
		      snprintf(par_name[numpar+i], 250,
			       "lin_t_%s_%d", ser[s].name,i);
		    }
		}
	      
	      numpar+=2; // update the number of parameters
	    }
	  else // global linear time dependency?
	    {
	      if(pars) // parameter value array is given?
		// fetch the parameter value:
		for(i=0;i<numsites;i++) // traverse the sites
		  ser[s].lin_t[i]=pars[numpar];
	      if(!pars)
		{
		  // fill the global parameter arrays with
		  // appropriate contents:
		  par_trans_type[numpar]=T_LIN;
		  par_type[numpar]=LIN_T;
		  par_layer[numpar]=ser[s].numlayers+1;
		  par_region[numpar]=-1;
		  par_series[numpar]=s;
		  snprintf(par_name[numpar], 250, "lin_t_%s",ser[s].name);
		}
	      
	      numpar++; // update the number of parameters
	    }
	}
      

      if(!ser[s].no_layers) // positive number of layers?
	{
	  // Fetch layer specific parameters:
	  for(l=0;l<ser[s].numlayers;l++)
	    if(!ser[s].time_integral[l])
	      {
		// handle the pull of this layer:
		if(l==(ser[s].numlayers-1) && ser[s].no_pull_lower) // no pull?
		  {
		    // set the pull extremly low:
		    for(i=0;i<numsites;i++)
		      if(ser[s].init_treatment)
			ser[s].pull[l][i]=0.0;
		      else
			ser[s].pull[l][i]=1e-10;
		  }
		else if(ser[s].regional_pull[l]) // regional pull?
		  {
		    for(i=0;i<numsites;i++) // traverse the sites
		      if(pars) // parameter value array is given?
			// set the parameter value:
			{
			  if(!ser[s].allow_positive_pulls)
			    ser[s].pull[l][i]=1.0/exp(pars[i+numpar]);
			  else
			    ser[s].pull[l][i]=pars[i+numpar];
			}
		      else
			{
			  // fill the global parameter arrays with
			  // appropriate contents:
			  if(!ser[s].allow_positive_pulls)
			    par_trans_type[i+numpar]=T_LOG_INV;
			  else
			    par_trans_type[i+numpar]=T_INV;
			  par_type[i+numpar]=DT;
			  par_layer[i+numpar]=l+1;
			  par_region[i+numpar]=i;
			  par_series[i+numpar]=s;
			  snprintf(par_name[i+numpar], 250,
				   "dt_%s_%d_%d",ser[s].name,l+1,i);
			}
		    
		    numpar+=numsites; // update the number of parameters
		  }
		else if(ser[s].indicator_pull[l]) // indicator-separated pull?
		  {
		    if(pars) // parameter value array is given?
		      {
			// set the parameter value:
			for(i=0;i<numsites;i++) // traverse the sites
			  if(!ser[s].allow_positive_pulls)
			    {
			      if(!indicator_array[i])
				ser[s].pull[l][i]=1.0/exp(pars[numpar]);
			      else
				ser[s].pull[l][i]=1.0/exp(pars[numpar+1]);
			    }
			  else
			    {
			      if(!indicator_array[i])
				ser[s].pull[l][i]=pars[numpar];
			      else
				ser[s].pull[l][i]=pars[numpar+1];
			    }
		      }
		    else
		      {
			for(i=0;i<2;i++) // traverse the indicator values
			  {
			    // fill the global parameter arrays with
			    // appropriate contents:
			    if(!ser[s].allow_positive_pulls)
			      par_trans_type[numpar+i]=T_LOG_INV;
			    else
			      par_trans_type[numpar+i]=T_INV;
			    par_type[numpar+i]=DT;
			    par_layer[numpar+i]=l+1;
			    par_region[numpar+i]=-1;
			    par_series[numpar+i]=s;
			    snprintf(par_name[numpar+i],250, 
				     "dt_%s_%d_I%d",ser[s].name, l+1,i);
			  }
		      }
		    
		    numpar+=2;  // update the number of parameters
		  }
		else // global pull
		  {
		    for(i=0;i<numsites;i++) // traverse the sites
		      if(pars) // parameter value array is given?
			// set the parameter value:
			{
			  if(!ser[s].allow_positive_pulls)
			    ser[s].pull[l][i]=1.0/exp(pars[numpar]);
			  else
			    ser[s].pull[l][i]=pars[numpar];
			}
		    if(!pars) // parameter value array is not given?
		      {
			// fill the global parameter arrays with
			// appropriate contents:
			if(!ser[s].allow_positive_pulls)
			  par_trans_type[numpar]=T_LOG_INV;
			else
			  par_trans_type[numpar]=T_INV;
			par_type[numpar]=DT;
			par_layer[numpar]=l+1;
			par_region[numpar]=-1;
			par_series[numpar]=s;
			snprintf(par_name[numpar], 250,
				 "dt_%s_%d",ser[s].name,l+1);
		      }
		    
		    numpar++; // update the number of parameters
		  }
		
		if(!ser[s].no_sigma[l]) // Diffusion on this layer?
		  // (stochastic layer?)
		  {
		    if(ser[s].regional_sigma[l]) // regional diffusion?
		      {
			for(i=0;i<numsites;i++) // traverse the sites
			  if(pars) // parameter value array is given?
			    // set the parameter value:
			    ser[s].sigma[l][i]=exp(pars[i+numpar]);
			  else
			    {
			      // fill the global parameter arrays with
			      // appropriate contents:
			      par_trans_type[i+numpar]=T_LOG;
			      par_type[i+numpar]=SIGMA;
			      par_layer[i+numpar]=l+1;
			      par_region[i+numpar]=i;
			      par_series[i+numpar]=s;
			      snprintf(par_name[i+numpar], 250,
				       "sigma_%s_%d_%d",ser[s].name,l+1,i);
			    }
			
			numpar+=numsites; // update the number of parameters
		      }
		    else if(ser[s].indicator_sigma[l]) // indicator-separated diffusion
		      {
			if(pars) // parameter value array is given?
			  {
			    // set the parameter value:
			    for(i=0;i<numsites;i++) // traverse the sites
			      if(!indicator_array[i])
				ser[s].sigma[l][i]=exp(pars[numpar]);
			      else
				ser[s].sigma[l][i]=exp(pars[numpar+1]);
			  }
			else  
			  {
			    for(i=0;i<2;i++) // traverse the indicator values
			      {
				// fill the global parameter arrays with
				// appropriate contents:
				par_trans_type[numpar+i]=T_LOG;
				par_type[numpar+i]=SIGMA;
				par_layer[numpar+i]=l+1;
				par_region[numpar+i]=-1;
				par_series[numpar+i]=s;
				snprintf(par_name[numpar+i], 250,
					 "sigma%s_%d_I%d",
					 ser[s].name, l+1,i);
			      }
			  }
			numpar+=2; // update the number of parameters
		      }
		    else // globally defined diffusion
		      {
			for(i=0;i<numsites;i++) // traverse the sites
			  if(pars) // parameter value array is given?
			    // set the parameter value:
			    ser[s].sigma[l][i]=exp(pars[numpar]);
			if(!pars)
			  {
			    // fill the global parameter arrays with
			    // appropriate contents:
			    par_trans_type[numpar]=T_LOG;
			    par_type[numpar]=SIGMA;
			    par_layer[numpar]=l+1;
			    par_region[numpar]=-1;
			    par_series[numpar]=s;
			    snprintf(par_name[numpar], 250,
				     "sigma_%s_%d",
				    ser[s].name,l+1);
			  }
			
			numpar++; // update the number of parameters
		      }
		  }
		else // no diffusion (deterministic upper layer)
		  for(i=0;i<numsites;i++)
		    {
		      if(!do_realize)
			ser[s].sigma[l][i]=0.0;
		      else
			{
			  ser[s].sigma[l][i]=ser[s].pr->s_1;
			  //cout << "sigma0=" << ser[s].sigma[l][i] << endl;
			}
		    }
		
		if(ser[s].sigma_pairwise_correlated[l])
		  {
		    if(pars)  // parameter value array is given?
		      {
			// set the parameter values:
			for(i=0;i<numsites;i++)
			  ser[s].paircorr[l][i][i]=1.0;
			
			for(i=0;i<(numsites-1);i++)			
			  for(j=i+1;j<numsites;j++)
			    ser[s].paircorr[l][i][j]=ser[s].paircorr[l][j][i]=
			      invtransform_parameter(pars[numpar++],T_LOGIST_PAIR);


      
			
			double *sigma_lambda=double_eigenvalues(ser[s].paircorr[l],
								numsites);
	
      

			for(i=0;i<numsites;i++)
			  if(sigma_lambda[i]<0.0 || !(sigma_lambda[i]<1e+200))
			    {
			      /* R debug code
				 cout << "sigma=matrix(c(";
				 for(i=0;i<numsites;i++)
				 for(j=0;j<numsites;j++)
				 {
				 cout << ser[s].paircorr[l][i][j];
				 if(i<(numsites-1) || j<(numsites-1))
				 cout << " , ";
				 }
				 cout << "),nrow=" << numsites << ")" << endl;
				 
				 cout << "lambda=c(";
				 for(i=0;i<numsites;i++)
				 {
				 cout << sigma_lambda[i];
				 if(i<(numsites-1))
				 cout << ",";
				 }
				 cout << ")" << endl;
			      */
			      
			      pairwise_wrong=true;
			    }
			
			delete [] sigma_lambda;
		      }
		    else
		      {
			for(i=0;i<(numsites-1);i++)
			  for(j=i+1;j<numsites;j++)
			    {
			      // fill the global parameter arrays with
			      // appropriate contents:
			      par_trans_type[numpar]=T_LOGIST_PAIR;
			      par_type[numpar]=PAIR_CORR;
			      par_layer[numpar]=l+1;
			      par_region[numpar]=i;
			      par_series[numpar]=s;
			      snprintf(par_name[numpar],250, 
				       "corr_%s_%d_%d,%d",
				       ser[s].name,l+1,i,j);
			      numpar++;
			    }
		      }
		  }
		
		if(ser[s].sigma_correlated[l])
		  {
		    if(pars)  // parameter value array is given?
		      // set the parameter value:
		      ser[s].corr[l]=invtransform_parameter(pars[numpar],T_LOGIST_GLOBAL);
		    else
		      {
			// fill the global parameter arrays with
			// appropriate contents:
			par_trans_type[numpar]=T_LOGIST_GLOBAL;
			par_type[numpar]=CORR;
			par_layer[numpar]=l+1;
			par_region[numpar]=-1;
			par_series[numpar]=s;
			snprintf(par_name[numpar], 250,
				 "corr_%s_%d",ser[s].name,l+1);
		      }
		    
		    numpar++; // update the number of parameters
		  }
		else // no inter-regional correlations in the upper layer?
		  // set the correlations to zero:
		  ser[s].corr[l]=0.0;
	      }
	    else // time integral
	      for(i=0;i<numsites;i++)
		ser[s].sigma[l][i]=ser[s].pull[l][i]=0.0;
	}
      else // no layers
	{
	  if(pars)
	    {
	      // Fix this by having one layer with hard pull and no diffusion
	      
	      for(i=0;i<numsites;i++)
		{
		  ser[s].pull[0][i]=0.001/ser[s].pr->dt_1;
		  ser[s].sigma[0][i]=0.0;
		}
	    }
	}
      
      int exist_missing_sd=0;
      for(k=0;k<ser[s].meas_len;k++)
	if(ser[s].meas[k].sd==MISSING_VALUE)
	  exist_missing_sd=1;
      
      if(exist_missing_sd)
	{
	  if(ser[s].pr->los_m==MISSING_VALUE || 
	     ser[s].pr->los_s==MISSING_VALUE)
	    {
#ifdef MAIN
	      cerr << "Observational noise prior needed but not "
		"given for serie \"" << ser[s].name << "\"!" << endl;
	      exit(0);
#else
	      Rcout << "Observational noise prior needed but not "
		"given for serie \"" << ser[s].name << "\"!" << std::endl;
#endif // MAIN
	    }
	  
	  if(pars)
	    ser[s].obs_sd=exp(pars[numpar]);
	  else
	    {
	      // fill the global parameter arrays with
	      // appropriate contents:
	      par_trans_type[numpar]=T_LOG;
	      par_type[numpar]=OBS_SD;
	      par_layer[numpar]=0;
	      par_region[numpar]=-1;
	      par_series[numpar]=s;
	      snprintf(par_name[numpar], 250, "obs_sd_%s",ser[s].name);
	    }
	  
	  numpar++;
	}
      
      if(ser[s].num_per>0)
	{
	  for(i=0;i<ser[s].num_per;i++) // traverse the periods
	    {
	      if(pars) // parameter value array is given?
		// fetch the parameter values:
		{
		  ser[s].beta_sin[i]=pars[numpar];
		  numpar++;
		  ser[s].beta_cos[i]=pars[numpar];
		  numpar++;
		}
	      else
		{
		  // Set global parameter arrays with
		  // appropriate contents:
		  par_trans_type[numpar]=T_LIN;
		  par_type[numpar]=TRIG;
		  par_layer[numpar]=ser[s].numlayers+1;
		  par_region[numpar]=-1;
		  par_series[numpar]=s;
		  snprintf(par_name[numpar], 250,
			   "beta_sin%d_%s",i+1,ser[s].name);
		  numpar++; // update the number of parameters
		  // Set global parameter arrays with
		  // appropriate contents:
		  
		  par_trans_type[numpar]=T_LIN;
		  par_type[numpar]=TRIG;
		  par_layer[numpar]=ser[s].numlayers+1;
		  par_region[numpar]=-1;
		  par_series[numpar]=s;
		  snprintf(par_name[numpar], 250,
			   "beta_cos%d_%s",i+1,ser[s].name);
		  numpar++; // update the number of parameters
		}
	    }
	}

    }

#ifndef MAIN
  if(detailed_loglik)
    Rcout << "start useext" << std::endl;
#endif //MAIN
  
  if(useext) // external timeseries 
    {
      if(pars) // parameter value array is given?
	// set the parameter value:
	ser[0].beta=pars[numpar];
      else
	{
	  // fill the global parameter arrays with
	  // appropriate contents:
	  par_trans_type[numpar]=T_LIN;
	  par_type[numpar]=BETA;
	  par_layer[numpar]=2;
	  par_region[numpar]=-1;
	  par_series[numpar]=0;
	  snprintf(par_name[numpar], 99, "beta");
	}
      
      numpar++; // update the number of parameters
    }
  else
    ser[0].beta=0.0; // set the regression coefficient to zero

#ifndef MAIN
  if(detailed_loglik)
    Rcout << "end useext, start init treatment" << std::endl;
#endif // MAIN
  
  for(s=0;s<num_series;s++)
    if(!ser[s].no_layers)
      {
	if(ser[s].init_treatment && ser[s].init_value==MISSING_VALUE)
	  {
	    if(ser[s].regional_init && ser[s].layered_init)
	      {
		for(l=0;l<ser[s].numlayers;l++)
		  for(i=0;i<numsites;i++)
		    {
		      if(pars)
			ser[s].init_mu[l][i]=pars[numpar];
		      else
			{
			  par_trans_type[numpar]=T_LIN;
			  par_type[numpar]=INIT;
			  par_layer[numpar]=l;
			  par_region[numpar]=i;
			  par_series[numpar]=s;
			  snprintf(par_name[numpar], 250,
				   "init_%s_l%d_s%d",
				   ser[s].name,l+1,i);
			}
		      numpar++;
		    }
	      }
	    else if(ser[s].regional_init && !ser[s].layered_init)
	      {
		for(i=0;i<numsites;i++)
		  {
		    if(pars)
		      {
			for(l=0;l<ser[s].numlayers;l++)
			  ser[s].init_mu[l][i]=pars[numpar];
		      }
		    else
		      {
			par_trans_type[numpar]=T_LIN;
			par_type[numpar]=INIT;
			par_layer[numpar]=0;
			par_region[numpar]=i;
			par_series[numpar]=s;
			snprintf(par_name[numpar], 250,
				 "init_%s_s%d",
				 ser[s].name,i);
		      }
		    numpar++;
		  }
	      }
	    else if(!ser[s].regional_init && ser[s].layered_init)
	      {
		for(l=0;l<ser[s].numlayers;l++)
		  {
		    if(pars)
		      {
			for(i=0;i<numsites;i++)
			  ser[s].init_mu[l][i]=pars[numpar];
		      }
		    else
		      {
			par_trans_type[numpar]=T_LIN;
			par_type[numpar]=INIT;
			par_layer[numpar]=l;
			par_region[numpar]=-1;
			par_series[numpar]=s;
			snprintf(par_name[numpar], 250,
				 "init_%s_l%d",
				 ser[s].name,l+1);
		      }
		    numpar++;
		  }
	      }
	    else // !ser[s].regional_init && !ser[s].layered_init
	      {
		if(pars)
		  {
		    for(l=0;l<ser[s].numlayers;l++)
		      for(i=0;i<numsites;i++)
			ser[s].init_mu[l][i]=pars[numpar];
		  }
		else
		  {
		    par_trans_type[numpar]=T_LIN;
		    par_type[numpar]=INIT;
		    par_layer[numpar]=0;
		    par_region[numpar]=-1;
		    par_series[numpar]=s;
		    snprintf(par_name[numpar], 250,
			     "init_%s",ser[s].name);
		  }
		numpar++;
	      }
	  }
	else if(ser[s].init_treatment && ser[s].init_value!=MISSING_VALUE)
	  {
	    for(l=0;l<ser[s].numlayers;l++)
	      for(i=0;i<numsites;i++)
		ser[s].init_mu[l][i]=ser[s].pr->init_m;
	  }
      }

#ifndef MAIN
  if(detailed_loglik)
    Rcout << "start correlations" << std::endl;
#endif // MAIN
  
  for(i=0;i<num_series_corr;i++)
    {
      if(pars)  // parameter value array is given?
	// set the parameter value:
	series_corr[i]=invtransform_parameter(pars[numpar],T_LOGIST_GLOBAL);
      else
	{
	  // fill the global parameter arrays with
	  // appropriate contents:
	  par_trans_type[numpar]=T_LOGIST_PAIR;
	  par_type[numpar]=SERIES_CORR;
	  par_layer[numpar]=corr_from_layer[i];
	  par_region[numpar]=-1;
	  par_series[numpar]=corr_from_series[i];
	  snprintf(par_name[numpar], 250, "corr_%s,%d_%s,%d",
		  ser[corr_from_series[i]].name,corr_from_layer[i]+1,
		  ser[corr_to_series[i]].name,corr_to_layer[i]+1);
	}
      
      numpar++; // update the number of parameters
    }

#ifndef MAIN
  if(detailed_loglik)
    Rcout << "start feedback" << std::endl;
#endif // MAIN
  
  for(i=0;i<num_series_feed;i++)
    {
      if(pars)  // parameter value array is given?
	// set the parameter value:
	beta_feed[i]=invtransform_parameter(pars[numpar],T_LIN);
      else
	{
	  // fill the global parameter arrays with
	  // appropriate contents:
	  par_trans_type[numpar]=T_LIN;
	  par_type[numpar]=BETA;
	  par_layer[numpar]=feed_to_layer[i];
	  par_region[numpar]=-1;
	  par_series[numpar]=feed_to_series[i];
	  snprintf(par_name[numpar], 250, "beta_%s,%d_%s,%d",
		  ser[feed_from_series[i]].name,
		  feed_from_layer[i]+1,ser[feed_to_series[i]].name,
		  feed_to_layer[i]+1);
	}
      
      numpar++; // update the number of parameters
    }

#ifndef MAIN
  if(detailed_loglik)
    Rcout << "timer stuff starting" << std::endl;
#endif // MAIN
  
  // If no parameter values were given, the global
  // parameter arrays should now have been filled and
  // it's time to return with the number of parameters 
  // as output:
  if(!pars) 
    {
#ifdef DETAILED_TIMERS
      timers[1][1]=clock();
      timers[1][2]+=(timers[1][1]-timers[1][0]);
      timers[2][1]=clock();
      timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
      return (double) numpar; 
    }
  
  
  
  // Sanity check on the parameters:
  if(pairwise_wrong)
    {
#ifdef DETAILED_TIMERS
      timers[1][1]=clock();
      timers[1][2]+=(timers[1][1]-timers[1][0]);
      timers[2][1]=clock();
      timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
      return -1e+200;
    }
  
  for(s=0;s<num_series;s++)
    for(i=0;i<numsites;i++)
      { 
	if(!(ser[s].mu[i]>-1e+200 && ser[s].mu[i]<1e+200))
	  {
#ifdef DETAILED_TIMERS
	    timers[1][1]=clock();
	    timers[1][2]+=(timers[1][1]-timers[1][0]);
	    timers[2][1]=clock();
	    timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
	    return -1e+200;
	  }
	for(l=0;l<ser[s].numlayers;l++)
	  if(!(ser[s].sigma[l][i]>-1e+200 && ser[s].sigma[l][i]<1e+200))
	    {
#ifdef DETAILED_TIMERS
	      timers[1][1]=clock();
	      timers[1][2]+=(timers[1][1]-timers[1][0]);
	      timers[2][1]=clock();
	      timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
	      return -1e+200;
	    }
	for(l=0;l<ser[s].numlayers;l++)
	  if(!(ser[s].pull[l][i]>-1e+200 && ser[s].pull[l][i]<1e+200))
	    {
#ifdef DETAILED_TIMERS
	      timers[1][1]=clock();
	      timers[1][2]+=(timers[1][1]-timers[1][0]);
	      timers[2][1]=clock();
	      timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
	      return -1e+200;
	    }
	for(l=0;l<ser[s].numlayers;l++)
	  if(!(ser[s].corr[l]>-1e+200 && ser[s].corr[l]<1e+200))
	    {
#ifdef DETAILED_TIMERS
	      timers[1][1]=clock();
	      timers[1][2]+=(timers[1][1]-timers[1][0]);
	      timers[2][1]=clock();
	      timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
	      return -1e+200;
	    }
	
	if(id_strategy!=ID_NONE)
	  // Is a layer tracking a process slower than the speed 
	  // of the process it is tracking? Do not accept this!
	  for(l=0;l<ser[s].numlayers-1;l++)
	    {
	      for(i=0;i<numsites;i++)
		if(ser[s].pull[l][i]<ser[s].pull[l+1][i] &&
		   !ser[s].time_integral[l])
		  {
#ifdef DETAILED_TIMERS
		    timers[1][1]=clock();
		    timers[1][2]+=(timers[1][1]-timers[1][0]);
		    timers[2][1]=clock();
		    timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
		    return -1e+200;
		  }
	    }
      }
  if(!(ser[0].beta>-1e+200 && ser[0].beta<1e+200))
    {
#ifdef DETAILED_TIMERS
      timers[1][1]=clock();
      timers[1][2]+=(timers[1][1]-timers[1][0]);
      timers[2][1]=clock();
      timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
      return -1e+200;
    }

#ifndef MAIN
  if(detailed_loglik)
    Rcout << "correlation2 starting" << std::endl;
#endif // MAIN
  
  //Rcout << "loglik corrtest" << std::endl;
  if(num_series_corr>0)
    {
      // Test if explicite series+layer correlations are 
      // positively definite
      
      double **sigma_series=Make_matrix(num_tot_layers,num_tot_layers);
      for(i=0;i<num_tot_layers;i++)
	sigma_series[i][i]=1.0;

      for(i=0;i<num_series_corr;i++)
	sigma_series[corr_from_index[i]][corr_to_index[i]]=
	  sigma_series[corr_to_index[i]][corr_from_index[i]]=series_corr[i];

      if(!check_matrix(sigma_series, num_tot_layers, num_tot_layers))
	{
#ifdef DETAILED_TIMERS
	  timers[1][1]=clock();
	  timers[1][2]+=(timers[1][1]-timers[1][0]);
	  timers[2][1]=clock();
	  timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
	  return -1e+200;

	}
      
      double *sigma_lambda=double_eigenvalues(sigma_series,num_tot_layers);
  		
      bool series_wrong=false;
      for(i=0;i<num_tot_layers;i++)
	if(sigma_lambda[i]<0.0 || !(sigma_lambda[i]<1e+200))
	  series_wrong=true;		    
			
      delete [] sigma_lambda;

      if(series_wrong)
	{
#ifdef DETAILED_TIMERS
	  timers[1][1]=clock();
	  timers[1][2]+=(timers[1][1]-timers[1][0]);
	  timers[2][1]=clock();
	  timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
	  return -1e+200;
	}
    }


  if(nodata)
    {
#ifdef DETAILED_TIMERS
      timers[1][1]=clock();
      timers[1][2]+=(timers[1][1]-timers[1][0]);
      timers[2][1]=clock();
      timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
      return 0.0;
    }
  
#ifdef MAIN
  FILE **f=NULL;
  if(simulation_files)
    {
      f=new FILE*[num_series];
      for(s=0;s<num_series;s++)
	f[s]=fopen(simulation_files[s],"w");
    }
#endif // MAIN
  

  
  // If we have arrived at this point, a likelihood is to be
  // calculated, using the Kalman filter approach...
  // (... or a simulation is to be performed)
  
  // Matrices and vectors used in the Kalman filter
  // (assumes numsites sites and numlayers=sum_series 
  // series.numlayers*numsites state variables):

#ifndef MAIN
  if(detailed_loglik)
    Rcout << "filling in matrices" << std::endl;
#endif // MAIN
  
  double **var=Make_matrix(num_states,num_states); // diffusion matrix

  double **A=Make_matrix(num_states,num_states);
  for(i=0;i<num_states;i++)
    {
      s=state_series[i];
      l=state_layer[i];
      
      if(!ser[s].time_integral[l])
	{
	  A[i][i] = -ser[s].pull[l][i%numsites];
	  if(l<(ser[s].numlayers-1))
	    A[i][i+numsites] = ser[s].pull[l][i%numsites];
	}
      else // time integral
	A[i][i+numsites] = 1.0;
    }
  
  is_complex=false;
  if(num_series_feed>0)
    {
      for(i=0;i<num_series_feed;i++)
	{
	  int i1=series_state_start[feed_from_series[i]]+numsites*feed_from_layer[i];
	  int i2=series_state_start[feed_to_series[i]]+numsites*feed_to_layer[i];   
	      
	  for(j=0;j<numsites;j++)
	    {
	      if(A[i2+j][i1+j]!=0.0)
		{
#ifdef MAIN
		  printf("Causal link from %s, layer %d to %s, layer %d found twice!",
			 ser[feed_from_series[i]].name,feed_from_layer[i]+1,
			 ser[feed_to_series[i]].name,feed_to_layer[i]+1);
		  exit(0);
#else
		  Rcout << "Causal link from " << ser[feed_from_series[i]].name <<
		    ", layer " << feed_from_layer[i]+1 << " to " <<
		    ser[feed_to_series[i]].name << ", layer " <<
		    feed_to_layer[i]+1 << " found twice!" << std::endl;
#endif // MAIN
		}
	      
	      A[i2+j][i1+j]=
		beta_feed[i]*ser[feed_to_series[i]].pull[feed_to_layer[i]][j];
	      
	      if(feed_symmetric[i])
		A[i1+j][i2+j]=A[i2+j][i1+j];
	    }
	}
    }
  
#ifndef MAIN
  if(detailed_loglik)
    Rcout << "eigenvector stuff" << std::endl;
#endif // MAIN

  // V=eigenvector matrix, Vinv=V^{-1},
  Complex *lambda=NULL, **V=NULL, **Vinv=NULL, **VinvLambdaV_A=NULL;
  double *lambda_r=NULL, **V_r=NULL, **Vinv_r=NULL;
  int numsites2=numsites; // for traversal of matrixes

  int allow_positive_pulls=ser[0].allow_positive_pulls;
  for(s=1;s<num_series;s++)
    allow_positive_pulls=allow_positive_pulls && ser[s].allow_positive_pulls;
  
  //Rcout << "loglik init" << std::endl;
  if(num_series_feed>0)
    {
      //V=get_complex_eigenvector_matrix(A,num_states,&lambda,0,1000000);
      //Vinv=inverse_complex_matrix(V,num_states);

      if(!check_matrix(A,num_states,num_states))
	{
#ifdef DETAILED_TIMERS
	  timers[1][1]=clock();
	  timers[1][2]+=(timers[1][1]-timers[1][0]);
	  timers[2][1]=clock();
	  timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
	  return -1e+200;
	}
      Vinv=matrix_eigenvectors(A,num_states,&lambda,&V);
      
      for(j=0;j<10;j++)
	cycle_time[j]=MISSING_VALUE;
      
      j=0;
      double prev_cycle=MISSING_VALUE;
      for(i=0;i<num_states && j<10 && j<ser[0].numlayers;i++)
	{
	  double cycle=MISSING_VALUE;
	  if(lambda[i].Im()!=0.0)
	    cycle=2*M_PI/ABSVAL((lambda[i].Im()));
	  if(cycle!=MISSING_VALUE && cycle!=prev_cycle)
	    {
	      cycle_time[j]=cycle;
	      //cout << cycle << " ";
	      j++;
	    }
	  prev_cycle=cycle;
	}
      //cout << " iscomplex=" << is_complex << " j=" << j << endl;
      
      for(i=0;i<num_states;i++)
	if(lambda[i].Re()>0.0 && !allow_positive_pulls)
	  {
	    //cout << "Positive pull:" << lambda[i].Re() << endl;
	    doubledelete(V,num_states);
	    doubledelete(Vinv,num_states);
	    delete [] lambda;
#ifdef DETAILED_TIMERS
	    timers[1][1]=clock();
	    timers[1][2]+=(timers[1][1]-timers[1][0]);
	    timers[2][1]=clock();
	    timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
	    return(-1e+200);
	  }
      
      // test the eigen-decomposition:
      VinvLambdaV_A=Make_Complex_matrix(num_states,num_states);
      for(i=0;i<num_states;i++)
	for(j=0;j<num_states;j++)
	  for(k=0;k<num_states;k++)
	    VinvLambdaV_A[i][j]+=Vinv[i][k]*lambda[k]*V[k][j];
      
      static int numtrav=0;
      for(i=0;i<num_states;i++)
	for(j=0;j<num_states;j++)
	  {
	    if(ABSVAL((VinvLambdaV_A[i][j].Re()-A[i][j]))>
	       (1e-3)*MAXIM(ABSVAL((VinvLambdaV_A[i][j].Re())),1.0) ||
	       ABSVAL((VinvLambdaV_A[i][j].Im()))>
	       (1e-3)*MAXIM(ABSVAL((VinvLambdaV_A[i][j].Re())),1.0))
	      {
		//cout << "VinvLambdaV_A[" << i << "," << j << "]=" << 
		//  VinvLambdaV_A[i][j] << "!=" << "A[" << i << "," << j << "]=" <<
		//  A[i][j] << endl;
		//cout << "counter=" << numtrav << endl;
#ifdef DETAILED_TIMERS
		timers[1][1]=clock();
		timers[1][2]+=(timers[1][1]-timers[1][0]);
		timers[2][1]=clock();
		timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
		return(-1e+200);
	      }
	  }
      numtrav++;
      
      for(i=0;i<num_states;i++)
	if(ABSVAL((lambda[i].Im()))>0.00001)
	  is_complex=true;
      if(!is_complex)
	{
	  for(i=0;i<num_states;i++)
	    for(j=0;j<num_states;j++)
	      if(ABSVAL((Vinv[i][j].Im()))>0.00001 ||
		 ABSVAL((V[i][j].Im()))>0.00001)
		is_complex=true;
	  
	  if(!is_complex)
	    {
	      lambda_r=new double[num_states];
	      Vinv_r=Make_matrix(num_states,num_states);
	      V_r=Make_matrix(num_states,num_states);

	      for(i=0;i<num_states;i++)
		{
		  lambda_r[i]=lambda[i].Re();
		  for(j=0;j<num_states;j++)
		    {
		      Vinv_r[i][j]=Vinv[i][j].Re();
		      V_r[i][j]=V[i][j].Re();
		    }
		}
	      
	      delete [] lambda;
	      lambda=NULL;
	      doubledelete(Vinv,num_states);
	      doubledelete(V,num_states);
	      doubledelete(VinvLambdaV_A,num_states);
	    }
	}
    }
  else
    {
      //V=get_complex_eigenvector_matrix(A,num_states,&lambda,0,1000000);
      //Vinv=inverse_complex_matrix(V,num_states);

      if(!check_matrix(A,num_states,num_states))
	{
#ifdef DETAILED_TIMERS
	  timers[1][1]=clock();
	  timers[1][2]+=(timers[1][1]-timers[1][0]);
	  timers[2][1]=clock();
	  timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
	  return -1e+200;
	}
      
      Vinv=matrix_eigenvectors(A,num_states,&lambda,&V);
      Vinv_r=Make_matrix(num_states,num_states);
      V_r=Make_matrix(num_states,num_states);
      lambda_r=new double[num_states];
      
      for(i=0;i<num_states;i++)
	{
	  lambda_r[i]=lambda[i].Re();
	  for(j=0;j<num_states;j++)
	    {
	      V_r[i][j]=V[i][j].Re();
	      Vinv_r[i][j]=Vinv[i][j].Re();
	    }
	}
      doubledelete(Vinv,num_states);
      doubledelete(V,num_states);
      delete [] lambda;
      lambda=NULL;

      for(i=0;i<num_states;i++)
	if(lambda_r[i]>0.0 && !allow_positive_pulls)
	  {
#ifdef MAIN
	    cout << "Positive pull:" << lambda_r[i] << endl;
#else
	    Rcout << "Positive pull:" << lambda_r[i] << std::endl;
#endif // MAIN
	    doubledelete(V_r,num_states);
	    doubledelete(Vinv_r,num_states);
	    delete [] lambda_r;
#ifdef DETAILED_TIMERS
	    timers[1][1]=clock();
	    timers[1][2]+=(timers[1][1]-timers[1][0]);
	    timers[2][1]=clock();
	    timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
  
	    return(-1e+200);
	  } 
    }

#ifndef MAIN
  if(detailed_loglik)
    Rcout << "m vector stuff" << std::endl;
#endif // MAIN
  
  double *m=new double[num_states]; // expectation vector
  for(i=0;i<num_states;i++)
    {
      s=state_series[i];
      l=state_layer[i];
      
      if(l<(ser[s].numlayers-1))
	m[i]=0.0;
      else
	m[i]=ser[s].mu[i%numsites];
      
      for(j=0;j<num_series_feed;j++)
	if(feed_to_series[j]==(int)s && feed_to_layer[j]==(int)l)
	  m[i] -= beta_feed[j]*ser[feed_from_series[j]].mu[i%numsites];
      
      m[i]*=ser[s].pull[l][i%numsites];
    }
  
  unsigned int len = dosmooth ? meas_smooth_len : meas_tot_len;
  
  measurement_cluster *me = (dosmooth ? meas_smooth : meas_tot);
  
#ifndef MAIN
  if(detailed_loglik)
    {
      Rcout << "dosmooth=" << dosmooth << std::endl;
      Rcout << "me=" << me << " meas_smooth=" << meas_smooth <<
	" meas_tot=" << meas_tot << std::endl;
      Rcout << " me[0].tm=" << me[0].tm << std::endl;
      Rcout << "omega/Lambda matrix stuff" << std::endl;
    }
#endif // MAIN
  
  // Omega=V*sigma^2*V', Vvar=V*var
  Complex **Omega=NULL, 
    **Vvar=NULL, 
    // Qbuffer=Vinv*Lambda_k:
    **Qbuffer=NULL, **Q_k_buff=NULL;
  // lambda_k=eLambda_k*omega/(eigenvalues), eLambda_k=e^(lambda*timediff):
  Complex **Lambda_k=NULL, 
    *eLambda_k=NULL;  
  // P_k_buffer=F_k*P_k_now:
  
  double **Omega_r=NULL,**Vvar_r=NULL,**Qbuffer_r=NULL; 
  double **Lambda_k_r=NULL,*eLambda_k_r=NULL;  
  
  if(is_complex)
    {
      Omega=Make_Complex_matrix(num_states,num_states);
      Vvar=Make_Complex_matrix(num_states,num_states);
      Qbuffer=Make_Complex_matrix(num_states,num_states); 
      Q_k_buff=Make_Complex_matrix(num_states,num_states); 
      Lambda_k=Make_Complex_matrix(num_states,num_states);
      eLambda_k=new Complex[num_states];
    }
  else
    {
      Omega_r=Make_matrix(num_states,num_states);
      Vvar_r=Make_matrix(num_states,num_states);
      Qbuffer_r=Make_matrix(num_states,num_states); 
      Lambda_k_r=Make_matrix(num_states,num_states);
      eLambda_k_r=new double[num_states];  
    }
  
#ifndef MAIN
  if(detailed_loglik)
    Rcout << "Making P, Q, x matrices" << std::endl;
#endif // MAIN
  
  // P_k_buffer=F_k*P_k_now:  
  double **P_k_buffer=Make_matrix(num_states,num_states);
  // Q_k = Qbuffer*Vinv=Vinv*Lambda_k*Vinv is the 
  // prediction uncertainty not due to uncertainty of the previous state
  double ***Q_k= new double**[len];
  for(k=0;k<len;k++)
    Q_k[k]=Make_matrix(num_states,num_states); 
  // x_k_prev=expectancy of the state given previous observations:
  double *x_k_prev=new double[num_states];
  
  double step=useext ? extdata[1].x-extdata[0].x : 0.0;
  
  
  // return value:
  double ret=0.0;
  
  double *t_k=new double[len];
  HydDateTime *dt_k=NULL;
  if(ser[0].meas[0].dt!=NoHydDateTime)
    dt_k=new HydDateTime[len];

  
#ifndef MAIN
  if(detailed_loglik)
    Rcout << "Making residual and prior expectancy matrices" << std::endl;
#endif // MAIN
  
  // count number of observed series:
  unsigned int numobs=num_series*numsites;
  double **resids=new double*[numobs];
  double **prior_expectancy=new double*[numobs];
  double *resids_time=new double[meas_tot_len];
  for(i=0;i<numobs;i++)
    {
      resids[i]=new double[len];
      prior_expectancy[i]=new double[len];
      for(k=0;k<len;k++)
	{
	  resids[i][k]=MISSING_VALUE;
	  prior_expectancy[i][k]=MISSING_VALUE;
	}
    }
  for(k=0;k<meas_tot_len;k++)
    resids_time[k]=meas_tot[k].tm;
  if(debug && return_residuals)
    for(k=0;k<meas_tot_len;k++)
      meas_tot[k].print();
  
  // smoothing information

#ifndef MAIN
  if(detailed_loglik)
    Rcout << "Making smoothing matrices" << std::endl;
#endif // MAIN
  
  // x_k_now=state expectancy given observations up until 
  // and including now: 
  double **x_k_now=new double*[len];
  // F_k=Vinv*eLambda_k*V is used for updating 
  // the state, x_k,(k-1) = F_k * x_(k-1),(k-1) + u_k
  double ***F_k=new double**[len];
  // u_k=the expectancy contributions to the state
  double **u_k=new double*[len];
  // P_k_now=variance of the state given observations until now:
  double ***P_k_now=new double**[len];
  // P_k_prev=variance of the state given previous observations:
  double ***P_k_prev=new double**[len];
  double ***P_k_s=NULL;
  double **x_k_s=NULL;
  double **C_k=NULL, **PAbuffer=NULL;


  if(dosmooth)
    {
      P_k_s=new double**[meas_smooth_len];
      x_k_s=new double*[meas_smooth_len];
      if(!t_k_smooth)
	t_k_smooth=new double[meas_smooth_len];
      if(!dt_k_smooth && ser[0].meas[0].dt!=NoHydDateTime)
	dt_k_smooth=new HydDateTime[meas_smooth_len];

      for(k=0;k<meas_smooth_len;k++)
	{
	  x_k_s[k]=new double[num_states];
	  P_k_s[k]=Make_matrix(num_states,num_states);
	}
      C_k=Make_matrix(num_states,num_states);
      PAbuffer=Make_matrix(num_states,num_states);
    }
  
  for(k=0;k<len;k++)
    {
      x_k_now[k]=new double[num_states];
      u_k[k]=new double[num_states];
      F_k[k]=new double*[num_states];
      P_k_now[k]=new double*[num_states];
      P_k_prev[k]=new double*[num_states];
      for(i=0;i<num_states;i++)
	{
	  F_k[k][i]=new double[num_states];
	  P_k_now[k][i]=new double[num_states];
	  P_k_prev[k][i]=new double[num_states];
	}
    }
  
#ifndef MAIN
  if(detailed_loglik)
    Rcout << "Correlation3 starting" << std::endl;
#endif // MAIN

  // Fill out correlation matrix:

  // Start with series+layer-wise correlation matrices:
  double ***corr_layer=new double**[num_tot_layers];
  double ***corr_sqrt_layer=NULL;
  if(num_series_corr>0)
    corr_sqrt_layer=new double**[num_tot_layers];
  
  for(i=0;i<num_tot_layers;i++)
    {
      corr_layer[i]=Make_matrix(numsites,numsites);
      if(num_series_corr>0)
	corr_sqrt_layer[i]=Make_matrix(numsites,numsites);
    }
  
  // Set one on the diagonals:
  for(i=0;i<num_tot_layers;i++)
    for(j=0;j<numsites;j++)
      corr_layer[i][j][j]=1.0;
  
  //Rcout << "loglik init2" << endl;
  
  // Fill out the off-diagonal terms of the
  // series+layer-wise correlation matrices:
  for(s=0;s<num_series;s++)
    for(l=0;l<ser[s].numlayers;l++)
      {
	i=(series_state_start[s]+numsites*l)/numsites;
      
	if(ser[s].sigma_1dim[l])
	  {
	    for(j=0;j<(numsites-1);j++)
	      for(k=j+1;k<numsites;k++)
		corr_layer[i][j][k]=corr_layer[i][k][j]=
		  dosmooth ? 0.999999 : 1.0;
	  }
	else if(ser[s].sigma_pairwise_correlated[l])
	  {
	    for(j=0;j<(numsites-1);j++)
	      for(k=j+1;k<numsites;k++)
		if((!ser[s].indicator_corr[l] || 
		    (indicator_array[j]==indicator_array[k])) &&
		   (!ser[s].indicator_corr2[l] || 
		    (indicator_array[k]==0 && indicator_array[j]==0)))
		  corr_layer[i][j][k]=corr_layer[i][k][j]=
		    ser[s].paircorr[l][j][k];
	  }
	else if(ser[s].sigma_correlated[l])
	  {
	    for(j=0;j<(numsites-1);j++)
	      for(k=j+1;k<numsites;k++)
		if((!ser[s].indicator_corr[l] || 
		    (indicator_array[j]==indicator_array[k])) &&
		   (!ser[s].indicator_corr2[l] || 
		    (indicator_array[k]==0 && indicator_array[j]==0)))
		  corr_layer[i][j][k]=corr_layer[i][k][j]=ser[s].corr[l];
	  }

	if(num_series_corr>0) // Find square root matrices?
	  {
	    if(ser[s].sigma_1dim[l])
	      {
		for(j=0;j<numsites;j++)
		  for(k=0;k<numsites;k++)
		    corr_sqrt_layer[i][j][k]=1.0/sqrt(double(numsites));
	      }
	    else
	      {
		
		Complex *sigma_lambda, **sigma_V;
		
		if(!check_matrix(corr_layer[i], numsites, numsites))
		  {
#ifdef DETAILED_TIMERS
		    timers[1][1]=clock();
		    timers[1][2]+=(timers[1][1]-timers[1][0]);
		    timers[2][1]=clock();
		    timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
		  return -1e+200;
		  }
		
		Complex **sigma_Vinv=matrix_eigenvectors(corr_layer[i],
							 numsites,
							 &sigma_lambda,&sigma_V);
		
		/*
		  complex *sigma_lambda;
		  complex **sigma_V=get_complex_eigenvector_matrix(corr_layer[i],numsites,
		  &sigma_lambda,0,1000000);
		  complex **sigma_Vinv=inverse_complex_matrix(sigma_V,numsites);
		*/
		
		if(!sigma_V)
		  {
#ifdef MAIN
		    cerr << "Failed to find eigenvalue decomposition of "
		      "correlation matrix for series " << s << " layer " << l << endl;
		    exit(0);
#else
		    Rcerr << "Failed to find eigenvalue decomposition of "
		      "correlation matrix for series " << s << " layer " <<
		      l << std::endl;
#endif // MAIN
		  }
		
		Complex **sqrtbuffer=Make_Complex_matrix(numsites,numsites);
		
		for(j=0;j<numsites;j++)
		  for(k=0;k<numsites;k++)
		    for(unsigned int k2=0;k2<numsites;k2++)
		      sqrtbuffer[j][k]+=
			sigma_Vinv[j][k2]*complex_sqrt(sigma_lambda[k2])*sigma_V[k2][k];
		
		for(j=0;j<numsites;j++)
		  for(k=0;k<numsites;k++)
		    {
		      if(ABSVAL((sqrtbuffer[j][k].Im()))>0.001)
			{
#ifdef MAIN
			  cerr << "Complex result when finding the square matrix of "
			    "the correlation matrix for series " << s << 
			    " layer " << l << endl;
			  exit(0);
#else
			  Rcout << "Complex result when finding the square "
			    "matrix of the correlation matrix for series " <<
			    s << " layer " << l << std::endl;
#endif // MAIN
			}

		      corr_sqrt_layer[i][j][k]=sqrtbuffer[j][k].Re();
		    }
		
		doubledelete(sqrtbuffer,numsites);
		doubledelete(sigma_V,numsites);
		doubledelete(sigma_Vinv,numsites);
		delete [] sigma_lambda;
	      }
	  }
      }
  
  // Fill out total correlation matrix:
  double **corr=Make_matrix(num_states,num_states);
  for(s=0;s<num_series;s++)
    for(l=0;l<ser[s].numlayers;l++)
      {
	i=(series_state_start[s]+numsites*l)/numsites;
	
	for(j=0;j<numsites;j++)
	  for(k=0;k<numsites;k++)
	    corr[i*numsites+j][i*numsites+k]=corr_layer[i][j][k];
      }
  //Rcout << "loglik init3" << std::endl;
  if(num_series_corr>0)
    {
      for(i=0;i<num_series_corr;i++)
	{
	  int index1=corr_from_index[i];
	  int index2=corr_to_index[i];
	  
	  double **sqrtprod=Make_matrix(numsites,numsites);
	  for(j=0;j<numsites;j++)
	    for(k=0;k<numsites;k++)
	      for(unsigned int k2=0;k2<numsites;k2++)
		sqrtprod[j][k]+=
		  corr_sqrt_layer[index1][j][k2]*corr_sqrt_layer[index2][k2][k];
	  
	  for(j=0;j<numsites;j++)
	    for(k=0;k<numsites;k++)
	      corr[index1*numsites+j][index2*numsites+k]=
		corr[index2*numsites+k][index1*numsites+j]=
		series_corr[i]*sqrtprod[j][k];
	  
	  doubledelete(sqrtprod,numsites);
	}

      if(!check_matrix(corr, num_states, num_states))
	{
#ifdef DETAILED_TIMERS
	  timers[1][1]=clock();
	  timers[1][2]+=(timers[1][1]-timers[1][0]);
	  timers[2][1]=clock();
	  timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
	  return -1+200;
	}
      Complex *sigma_lambda=Complex_eigenvalues(corr, num_states);
			
      bool corr_wrong=false;
      for(i=0;i<num_states;i++)
	if(sigma_lambda[i].Re()<0.0 || !(sigma_lambda[i].Re()<1e+200))
	  corr_wrong=true;		    
      
      delete [] sigma_lambda;

      if(corr_wrong)
	{
#ifdef MAIN
	  cerr << "Series and site correlation matrices ok, yet total "
	    "correlation matrix is not!" << endl;
	  exit(0);
#else
	  Rcout << "Series and site correlation matrices ok, yet total "
	    "correlation matrix is not!" << std::endl;
#endif // MAIN
	}
    }
  
  // ********************************************
  // Fill out the covariance matrix, Sigma^2:

#ifndef MAIN
  if(detailed_loglik)
    Rcout << "Sigma matrix starting" << std::endl;
#endif // MAIN
  
  // First fill out the standard deviation vector:
  double *sd_vector=new double[num_states];
  for(i=0;i<num_states;i++)
    {
      s=state_series[i];
      l=state_layer[i];
      j=i%numsites;

      sd_vector[i]=ser[s].sigma[l][j];
    }
  
  // Fill out the diffusion matric based on the standard deviation vector and
  // the correlation matrix:
  for(i=0;i<num_states;i++)
    for(j=0;j<num_states;j++)
      var[i][j]=corr[i][j]*sd_vector[i]*sd_vector[j];

  // cleanup:
  for(i=0;i<num_tot_layers;i++)
    {
      doubledelete(corr_layer[i],numsites);
      if(num_series_corr>0)
	doubledelete(corr_sqrt_layer[i],numsites);
    }
  delete [] corr_layer;
  if(num_series_corr>0)
    delete [] corr_sqrt_layer;
  doubledelete(corr,num_states);
  delete [] sd_vector;
  // Done filling out the diffusion matrix


#ifndef MAIN
  if(detailed_loglik)
    Rcout << "Finding V*Sigma^2" << std::endl;
#endif // MAIN
  
  // ******************
  // Find V*Sigma^2
  for(i=0;i<num_states;i++)
    for(j=0;j<num_states;j++)
      {
	if(is_complex)
	  {
	    Vvar[i][j]=0.0;
	    for(l=i%numsites2;l<num_states;l+=numsites2)
	      Vvar[i][j]+=V[i][l]*var[l][j];
	  }
	else
	  {
	    int start_i=num_series_feed>0 ? i%numsites : i; 
	    Vvar_r[i][j]=0.0;
	    for(l=start_i;l<num_states;l+=numsites2)
	      Vvar_r[i][j]+=V_r[i][l]*var[l][j];
	  }
      }


#ifndef MAIN
  if(detailed_loglik)
    Rcout << "Finding V*Sigma^2*V'" << std::endl;
#endif // MAIN
  
  // *************************
  // Find Omega=V*sigma^2*V'
  for(i=0;i<num_states;i++)
    for(j=0;j<num_states;j++)
      {
	if(is_complex)
	  {
	    Omega[i][j]=0.0;
	    for(l=j%numsites2;l<num_states;l+=numsites2)
	      Omega[i][j]+=Vvar[i][l]*(V[j][l].conjugate());
	  }
	else
	  {
	    Omega_r[i][j]=0.0;
	    for(l=j%numsites2;l<num_states;l+=numsites2)
	      Omega_r[i][j]+=Vvar_r[i][l]*V_r[j][l];
	  }
      }
	


#ifndef MAIN
  if(detailed_loglik)
    {
      Rcout << "Initialize x,u,P,F" << std::endl;
      Rcout << "num_states=" << num_states << " num_series=" << num_series << std::endl;
      for(s=0;s<num_series;s++)
	Rcout << "mean time " << s << ":" << ser[s].mean_time << std::endl;
      Rcout << "is_complex=" << is_complex << std::endl;
    }
#endif // MAIN
  
  // ********************************************
  // Initialization of x_k_k, u_k, P_k_k and F_k:
  t_k[0]=me[0].tm;
  Complex *W=NULL,*u_k_buff=NULL;
  double *W_r=NULL,*u_k_buff_r=NULL;
  
  if(is_complex)
    {
      W=new Complex[num_states];
      u_k_buff=new Complex[num_states];
    }
  else
    {
#ifndef MAIN
      if(detailed_loglik)
	Rcout << "num_states=" << num_states << std::endl;
#endif // MAIN
      W_r=new double[num_states];
#ifndef MAIN
      if(detailed_loglik)
	Rcout << "W initialized" << std::endl;
#endif // MAIN
      u_k_buff_r=new double[num_states];
    }
  
#ifndef MAIN
  if(detailed_loglik)
    Rcout << "W,u initialized" << std::endl;
#endif // MAIN
  
  double meantime=0.0;
  for(s=0;s<num_series;s++)
    meantime+=ser[s].mean_time;
  meantime/=double(num_series);

#ifndef MAIN
  if(detailed_loglik)
    Rcout << "meantime filled" << std::endl;
#endif // MAIN
  
#ifdef MAIN
  if(simulation_files)
    {
      cout << "A:" << endl;
      for(i=0;i<num_states;i++)
	{
	  cout << i << ": ";
	  for(j=0;j<num_states;j++)
	    cout << A[i][j] << " ";
	  cout << endl;
	}
      cout << endl;
      
      cout << "V:" << endl;
      for(i=0;i<num_states;i++)
	{
	  cout << i << ": ";
	  for(j=0;j<num_states;j++)
	    if(is_complex)
	      cout << V[i][j] << " ";
	    else
	      cout << V_r[i][j] << " ";
	  cout << endl;
	}
      cout << endl;
      
      cout << "Vinv:" << endl;
      for(i=0;i<num_states;i++)
	{
	  cout << i << ": ";
	  for(j=0;j<num_states;j++)
	    if(is_complex)
	      cout << Vinv[i][j] << " ";
	    else
	      cout << Vinv_r[i][j] << " ";
	  cout << endl;
	}
      cout << endl;
      
      cout << "Lambda: ";
      for(i=0;i<num_states;i++)
	if(is_complex)
	  cout << lambda[i] << " ";
	else
	  cout << lambda_r[i] << " ";
      cout << endl;

      cout << "m: ";
      for(i=0;i<num_states;i++)
	cout << m[i] << " ";
      cout << endl;
    }
#endif // MAIN

#ifndef MAIN
  if(detailed_loglik)
    Rcout << "Fill W, initialize x, u" << std::endl;
#endif // MAIN
  
  for(i=0;i<num_states;i++)
    {
      s=state_series[i];
      
      x_k_now[0][i]=0.0;
      u_k[0][i]=0.0;
      
      if(is_complex)
	W[i]=u_k_buff[i]=0.0;
      else
	W_r[i]=u_k_buff_r[i]=0.0;
      
      l=state_layer[i];
      
      if(!ser[s].no_pull_lower)
	{
	  if((is_complex && lambda[i]!=0.0) ||
	     (!is_complex && lambda_r[i]!=0))
	    {
	      unsigned int b;
	      for(b=i%numsites2;b<num_states;b+=numsites2)
		{
		  if(is_complex)
		    W[i] -= 1.0/lambda[i]*V[i][b]*m[b];
		  else
		    W_r[i] -= 1.0/lambda_r[i]*V_r[i][b]*m[b];
		}
	      
	      if(ser[s].linear_time_dep)
		{
		  double t00=me[0].tm;
		  unsigned int nl=ser[s].numlayers;
		  b=series_state_start[s]+numsites*(nl-1)+i%numsites;
	      
	          if(is_complex)
		    W[i] -= ser[s].pull[nl-1][i%numsites]/lambda[i]*V[i][b]*
		      ser[s].lin_t[i%numsites]*
		      ((t00-meantime)+1.0/lambda[i]);
		  else
		    W_r[i] -= ser[s].pull[nl-1][i%numsites]/lambda_r[i]*V_r[i][b]*
		      ser[s].lin_t[i%numsites]*
		      ((t00-meantime)+1.0/lambda_r[i]);
		}
	    }
	  else
	    {
	      if(is_complex)
		W[i]=0.0; // Assume fixed initial value
	      else
		W_r[i]=0.0; // Assume fixed initial value
	    }
	}
    }
  
#ifndef MAIN
  if(detailed_loglik)
    Rcout << "Fill u" << std::endl;
#endif // MAIN
  
  for(i=0;i<num_states;i++) // traverse the states
    {
      if(is_complex)
	{
	  u_k_buff[i]=0.0;
	  for(j=0;j<num_states;j++)
	    if(Vinv[i][j]!=0.0)
	      u_k_buff[i]+=Vinv[i][j]*W[j];    
	  
	  if(ABSVAL((u_k_buff[i].Im()))>0.001)
	    {
#ifdef MAIN
	      cerr << "Complex expectancy! u_k[" << k << "][" << i << "]=" << 
		u_k_buff[i].Re() << "+" << u_k_buff[i].Im() << "i" << endl;
	      
	      cout << "A=matrix(c(";
	      for(j=0;j<num_states;j++)
		for(k=0;k<num_states;k++)
		  {
		    cout << A[j][k];
		    if(j<(num_states-1) || k<(num_states-1))
		      cout << ",";
		  }
	      cout << "),nrow=" << num_states << ")" << endl;
	      
	      cout << "lambda=c(";
	      for(j=0;j<num_states;j++)
		{
		  cout << lambda[j];
		  if(j<(num_states-1))
		    cout <<",";
		}
	      cout << ")" << endl;
	      
	      cout << "m=c(";
	      for(j=0;j<num_states;j++)
		{
		  cout << m[i];
		  if(j<(num_states-1))
		    cout <<",";
		}
	      cout << ")" << endl;
	      
	      cout << "V=matrix(c(";
	      for(j=0;j<num_states;j++)
		for(k=0;k<num_states;k++)
		  {
		    cout << V[j][k];
		    if(j<(num_states-1) || k<(num_states-1))
		      cout << ",";
		  }
	      cout << "),nrow=" << num_states << ")" << endl;
	      
	      cout << "Vinv=matrix(c(";
	      for(j=0;j<num_states;j++)
		for(k=0;k<num_states;k++)
		  {
		    cout << Vinv[j][k];
		    if(j<(num_states-1) || k<(num_states-1))
		      cout << ",";
		  }
	      cout << "),nrow=" << num_states << ")" << endl;
	      
	      cout << "VinvLambdaV=matrix(c(";
	      for(j=0;j<num_states;j++)
		for(k=0;k<num_states;k++)
		  {
		    cout << VinvLambdaV_A[j][k];
		    if(j<(num_states-1) || k<(num_states-1))
		      cout << ",";
		  }
	      cout << "),nrow=" << num_states << ")" << endl;
#else
	      Rcout << "Complex expectancy! u_k[" << k << "][" << i << "]=" << 
		u_k_buff[i].Re() << "+" << u_k_buff[i].Im() << "i" << std::endl;
	      
#endif // MAIN	      

	      doubledelete(prior_expectancy, numobs);
	      doubledelete(resids, numobs);
	      delete [] resids_time;
#ifdef DETAILED_TIMERS
	      timers[1][1]=clock();
	      timers[1][2]+=(timers[1][1]-timers[1][0]);
	      timers[2][1]=clock();
	      timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
	      return(-1e+200);
	      //exit(0);
	    }

	  u_k[0][i]=u_k_buff[i].Re();
	}
      else // no complex eigenvalues
	{
	  u_k_buff_r[i]=0.0;
	  int start_j = num_series_feed>0 ? i%numsites2 : i;
	  for(j=start_j;j<num_states;j+=numsites2)
	    if(Vinv_r[i][j]!=0.0)
	      u_k_buff_r[i]+=Vinv_r[i][j]*W_r[j];    
	  u_k[0][i]=u_k_buff_r[i];
	}
    }
  
  if(simulation_files)
    {
#ifdef MAIN
      cout << "W_0: ";
      for(i=0;i<num_states;i++)
	if(is_complex)
	  cout << W[i] << " ";
	else
	  cout << W_r[i] << " ";
      cout << endl;
      
      cout << "u_k_0: ";
      for(i=0;i<num_states;i++)
	cout << u_k[0][i] << " ";
      cout << endl;
#endif // MAIN
    }

  for(i=0;i<num_states;i++)
    for(j=0;j<num_states;j++)
      {
	P_k_now[0][i][j]=0.0;
	F_k[0][i][j]=0.0;
      }
  

#ifndef MAIN
  if(detailed_loglik)
    Rcout << "useext2 starting" << std::endl;
#endif // MAIN
  
  // ************************************************************
  // If there is an external timeseries,
  // do numerical integration up until the first measurement time:
  //double t0=useext ? extdata[0].x : 0.0; // starting time 
  // of the external timeseries
  double t1; //time closest to the first measurement:
  if(useext)
    {
      // find the external timeseries time closest to the first
      // measurement:
      for(t=0;t<ext_len && (extdata[t].x+0.5*step)<(me[0].tm);t++);
      t1=extdata[t].x;

      // Find the time difference between the starting time of
      // the timeseries and that closest to the first measurement:
      double dt2=t1;

      for(i=0;i<num_states;i++)
	{
	  if(is_complex)
	    W[i]=0.0;
	  else
	    W_r[i]=0.0;
	}
      
      // Traverse the time steps in the external timeseries
      // previous to the first measurement:
      for(t=0;t<ext_len && (extdata[t].x+0.5*step)<(me[0].tm);t++)
	{
	  double dt=t>0 ? (extdata[t].x-extdata[t-1].x) : step;
	  double u=extdata[t].x;

	  for(i=0;i<ser[0].numlayers*numsites;i++)
	    if(is_complex)
	      {
		Complex W_t=0.0;
		
		for(j=numsites;j<2*numsites;j++)
		  if(V[i][j]!=0.0)
		    {
		      if(j>=numsites && j<2*numsites)
			{
			  Complex ee=lambda[i]*(dt2-u);
			  W_t+=complex_exp(ee)*V[i][j]*
			    extdata[t].y*ser[0].beta;
			}
		    }
		
		W[i]+=W_t*dt;
	      }
	    else
	      {
		double W_t=0.0;
		
		for(j=numsites;j<2*numsites;j++)
		  if(V_r[i][j]!=0.0)
		    {
		      if(j>=numsites && j<2*numsites)
			{
			  double ee=lambda_r[i]*(dt2-u);
			  W_t+=exp(ee)*V_r[i][j]*extdata[t].y*ser[0].beta;
			}
		    }
		
		W_r[i]+=W_t*dt;
	      }
	}
      
      for(i=0;i<num_states;i++) // traverse the sites
	if(is_complex)
	  {
	    u_k_buff[i]=u_k[0][i];
	    for(j=0;j<num_states;j++)
	      if(Vinv[i][j]!=0.0)
		u_k_buff[i]+=Vinv[i][j]*W[j];
	    
	    
	    if(ABSVAL((u_k_buff[i].Im()))>0.001)
	      {
#ifdef MAIN
		cerr << "Complex expectancy! u_k[" << k << "][" << i << "]=" << 
		  u_k_buff[i].Re() << "+" << u_k_buff[i].Im() << "i" << endl;
#else
		Rcout << "Complex expectancy! u_k[" << k << "][" <<
		  i << "]=" << u_k_buff[i].Re() << "+" <<
		  u_k_buff[i].Im() << "i" << std::endl;
#endif // MAIN
		doubledelete(prior_expectancy, numobs);
		doubledelete(resids, numobs);
		delete [] resids_time;
#ifdef DETAILED_TIMERS
		timers[1][1]=clock();
		timers[1][2]+=(timers[1][1]-timers[1][0]);
		timers[2][1]=clock();
		timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
		return(-1e+200);
	      }
	    
	    u_k[0][i]=u_k_buff[i].Re();
	  }
	else
	  {
	    u_k_buff_r[i]=u_k[0][i];
	    for(j=0;j<num_states;j++)
	      if(Vinv_r[i][j]!=0.0)
		u_k_buff_r[i]+=Vinv_r[i][j]*W_r[j];
	    u_k[0][i]=u_k_buff_r[i];

	    if(simulation_files && debug)
	      {
#ifdef MAIN
		double upd=0.0;
		cout << " W[" << i << "]=" << W_r[i] << endl;
		for(j=0;j<num_states;j++)
		  {
		    cout << Vinv_r[i][j] << " ";
		    if(Vinv_r[i][j]!=0.0)
		      upd+=Vinv_r[i][j]*W_r[j];
		  }
		cout << endl;
		cout << " u_k_upd" << i << ":" << upd << endl;
#endif // MAIN
	      }
	  }
      
      t_0=t;
      //t0=t1;
    }
  

  if(simulation_files && debug)
    {
#ifdef MAIN
      cout << "u_k_2: ";
      for(i=0;i<num_states;i++)
	cout << u_k[0][i] << " ";
      cout << endl;
#else
      Rcout << "u_k_2: ";
      for(i=0;i<num_states;i++)
	Rcout << u_k[0][i] << " ";
      Rcout << std::endl;
#endif // MAIN
    }

  
#ifdef DETAILED_TIMERS
  timers[2][1]=clock();
  timers[2][2]+=(timers[2][1]-timers[2][0]);
#endif // DETAILED_TIMERS
  
  
#ifdef DETAILED_TIMERS
  timers[3][0]=clock();
#endif // DETAILED_TIMERS

  if(detailed_loglik)
    {
      show_mat("sigma^2",var,num_states,num_states);
      if(is_complex)
	show_Complex_mat("Omega", Omega,num_states,num_states);
      else
	show_mat("Omega_r", Omega_r,num_states,num_states);
      
      show_mat("A",A,num_states,num_states);
      if(!is_complex)
	show_vec("lambda_r",lambda_r,num_states);
      else
	show_Complex_vec("lambda",lambda,num_states);
      if(!is_complex)
	show_mat("V_r",V_r,num_states,num_states);
      else
	show_Complex_mat("V",V,num_states,num_states);
      
      if(is_complex)
	show_Complex_mat("VinvLambdaV_A",VinvLambdaV_A,num_states,num_states);
    }
  
  // *********************************
  // Traverse the measurements:
  for(k=0;k<len;k++)
    {
#ifndef MAIN
      if(detailed_loglik)
	Rcout << "k=" << k << std::endl;
#endif // MAIN
  
      // fetch the time of the current measurement:
      t_k[k]=me[k].tm;
      
      if(ser[0].meas[0].dt!=NoHydDateTime)
	dt_k[k]=me[k].dt;
      
      if(k==0) // first?
	{
	  if(is_complex)
	    {
	      // calculate the Lambda matrix
	      for(i=0;i<num_states;i++) // traverse the state vector
		for(j=0;j<num_states;j++) // traverse the state vector
		  if(lambda[i]!=0.0 || lambda[j]!=0.0)
		    Lambda_k[i][j]=-1.0/(lambda[i]+(lambda[j].conjugate()))*Omega[i][j];
		  else
		    Lambda_k[i][j]=0.0;
	    }
	  else
	    {
	      // calculate the Lambda matrix
	      for(i=0;i<num_states;i++) // traverse the state vector
		for(j=0;j<num_states;j++) // traverse the state vector
		  if(lambda_r[i]!=0.0 || lambda_r[j]!=0.0)
		    Lambda_k_r[i][j]=-1.0/(lambda_r[i]+lambda_r[j])*Omega_r[i][j];
		  else
		    Lambda_k_r[i][j]=0.0;
	    }
	}
      else // not first
	{
	  // Calculate u_k:
	  double dt=t_k[k]-t_k[k-1];
	  
	  for(i=0;i<num_states;i++) // traverse the sites
	    {
	      u_k[k][i]=0.0;
	      if(is_complex)
		W[i]=0.0;
	      else
		W_r[i]=0.0;
	      
	      s=state_series[i];
	      l=state_layer[i];
	      unsigned int nl=ser[s].numlayers;
	      unsigned int b=i+(nl-l-1)*numsites;
	      
	      if(state_layer[b]!=(nl-1))
		{
#ifdef MAIN
		  cerr << "state_layer=" << state_layer[b] <<
		    " numlayers-1=" << nl-1 << endl;
		  exit(0);
#else
		  Rcout << "state_layer=" << state_layer[b] <<
		    " numlayers-1=" << nl-1 << std::endl;
#endif // MAIN
		}
	      
	      if(!ser[s].no_pull_lower &&
		 ((is_complex && lambda && lambda[i]!=0.0) ||
		  (!is_complex && lambda_r && lambda_r[i]!=0.0)))
		{
		  if(is_complex)
		    {
		      Complex ee=lambda[i]*dt;
		      
		      for(b=i%numsites2;b<num_states;b+=numsites2)
			W[i] -= 1.0/lambda[i]*V[i][b]*m[b]*(1.0-complex_exp(ee));
		    }
		  else
		    {
		      double ee=lambda_r[i]*dt;
		      
		      for(b=i%numsites2;b<num_states;b+=numsites2)
			W_r[i] -= 1.0/lambda_r[i]*V_r[i][b]*m[b]*(1.0-exp(ee));
		    }
		  
		  if(ser[s].linear_time_dep)
		    {
		      b=series_state_start[s]+numsites*(nl-1)+i%numsites;
		      
		      if(is_complex)
			{
			  Complex ee=lambda[i]*dt;
			  W[i] -= ser[s].pull[nl-1][i%numsites]/
			    lambda[i]*V[i][b]*ser[s].lin_t[i%numsites]*
			    ((t_k[k]-meantime)-
			     (t_k[k-1]-meantime)*complex_exp(ee)+
			     1.0/lambda[i]*(1.0-complex_exp(ee)));
			}
		      else
			{
			  double ee=lambda_r[i]*dt;
			  W_r[i] -= ser[s].pull[nl-1][i%numsites]/
			    lambda_r[i]*V_r[i][b]*ser[s].lin_t[i%numsites]*
			    ((t_k[k]-meantime)-
			     (t_k[k-1]-meantime)*exp(ee) +
			     1.0/lambda_r[i]*(1.0-exp(ee)));
			}
		    }
		}
	      
	      // If an external timeseries is given:
	      if(useext && state_series[i]==0)
		{
		  // Update the numerical integral:

		  // find the external timeseries time closest to the 
		  // current measurement:
		  for(t=t_0;t<ext_len && (extdata[t].x+0.5*step)<(me[k].tm);t++);
		  t1=extdata[t].x;
	      
		  // Traverse the time steps in the external timeseries
		  // from the previously handled time to the current measurement:
		  for(t=t_0;t<ext_len && (extdata[t].x+0.5*step)<(me[k].tm);t++)
		    {
		      // time difference between currently handled time and start time:
		      double u=extdata[t].x;
		      double dt=extdata[t].x-extdata[t-1].x;

		      if(is_complex)
			{
			  Complex W_t=0.0;
			  Complex ee=lambda[i]*(t1-u);
			  
			  for(j=numsites;j<2*numsites;j++)
			    if(V[i][j]!=0.0)
			      W_t+=complex_exp(ee)*V[i][j]*extdata[t].y;
			  
			  W[i]+=ser[0].beta*W_t*dt;
			}
		      else
			{
			  double W_t=0.0;
			  double ee=lambda_r[i]*(t1-u);
			  
			  for(j=numsites;j<2*numsites;j++)
			    if(V_r[i][j]!=0.0)
			      W_t+=exp(ee)*V_r[i][j]*extdata[t].y;
			  
			  W_r[i]+=ser[0].beta*W_t*dt;
			}
		    }
		}
	    }
	  if(useext)
	    t_0=t;
	  
	  for(i=0;i<num_states;i++) // traverse the sites
	    if(is_complex)
	      {
		u_k_buff[i]=0.0;
		for(j=i%numsites2;j<num_states;j+=numsites2)
		  if(Vinv[i][j]!=0.0)
		    u_k_buff[i]+=Vinv[i][j]*W[j];
		
		if(ABSVAL((u_k_buff[i].Im()))>0.001)
		  {
#ifdef MAIN
		    cerr << "Complex expectancy! u_k[" << k << "][" <<
		      i << "]=" << u_k_buff[i].Re() << "+" <<
		      u_k_buff[i].Im() << "i" << endl;
#else
		    Rcout << "Complex expectancy! u_k[" << k << "][" <<
		      i << "]=" << u_k_buff[i].Re() << "+" <<
		      u_k_buff[i].Im() << "i" << std::endl;
#endif // MAIN
		    doubledelete(prior_expectancy, numobs);
		    doubledelete(resids, numobs);
		    delete [] resids_time;
#ifdef DETAILED_TIMERS
		    timers[1][1]=clock();
		    timers[1][2]+=(timers[1][1]-timers[1][0]);
		    timers[3][1]=clock();
		    timers[3][2]+=(timers[3][1]-timers[3][0]);
#endif // DETAILED_TIMERS
		    return(-1e+200);
		  }
		
		u_k[k][i]=u_k_buff[i].Re();
	      }
	    else
	      {
		u_k_buff_r[i]=0.0;
		for(j=i%numsites2;j<num_states;j+=numsites2)
		  if(Vinv_r[i][j]!=0.0)
		    u_k_buff_r[i]+=Vinv_r[i][j]*W_r[j];
		
		u_k[k][i]=u_k_buff_r[i];
	      }

	  
	  // find e^{lambda*timediff}
	  for(i=0;i<num_states;i++)
	    if(is_complex)
	      {
		Complex ee=lambda[i]*dt;
		eLambda_k[i]=complex_exp(ee);
	      }
	    else
	      eLambda_k_r[i]=exp(lambda_r[i]*dt);
	  
	  // calculate the Lambda_k matrix
	  for(i=0;i<num_states;i++)
	    for(j=0;j<num_states;j++)
	      if(is_complex)
		{
		  if(Omega[i][j]!=0.0)
		    {
		      if(lambda[i]!=0.0 || lambda[j]!=0.0)
			Lambda_k[i][j]=(eLambda_k[i]*(eLambda_k[j].conjugate())-1.0)/
			  (lambda[i]+(lambda[j].conjugate()))*Omega[i][j];
		      else
			Lambda_k[i][j]=Omega[i][j]*dt;
		    }
		}
	      else
		{
		  if(Omega_r[i][j]!=0.0)
		    {
		      if(lambda_r[i]!=0.0 || lambda_r[j]!=0.0)
			Lambda_k_r[i][j]=(eLambda_k_r[i]*eLambda_k_r[j]-1.0)/
			  (lambda_r[i]+lambda_r[j])*Omega_r[i][j];
		      else
			Lambda_k_r[i][j]=Omega_r[i][j]*dt;
		    }
		}
	  
	  if(is_complex)
	    {
	      Complex **F_k_buff=Make_Complex_matrix(num_states,num_states);
	      // Calculate F_k=Vinv*e^lambda(t(k)-t(k-1))*V :
	      // F_k=Vinv*eLambda_k*V:
	      for(i=0;i<num_states;i++)
		for(j=i%numsites2;j<num_states;j+=numsites2)
		  for(l=i%numsites2;l<num_states;l+=numsites2)
		    F_k_buff[i][j]+=Vinv[i][l]*V[l][j]*eLambda_k[l];
	      // PS: F_k will also be a sparse upper triagonal matrix, just as V
	      // *if* there are no feedback loops
	      
	      for(i=0;i<num_states;i++)
		for(j=0;j<num_states;j++)
		  {
		    if(ABSVAL((F_k_buff[i][j].Im()))>0.001)
		      {
#ifdef MAIN
			cerr << "Complex linear transformation! F_k[" <<
			  k << "][" << 
			  i << "][" << j << "]=" << 
			  F_k_buff[i][j].Re() << "+" << 
			  F_k_buff[i][j].Im() << "i" << endl;
#else
			Rcout << "Complex linear transformation! F_k[" <<
			  k << "][" << 
			  i << "][" << j << "]=" << 
			  F_k_buff[i][j].Re() << "+" << 
			  F_k_buff[i][j].Im() << "i" << std::endl;
#endif // MAIN
			doubledelete(prior_expectancy, numobs);
			doubledelete(resids, numobs);
			delete [] resids_time;
#ifdef DETAILED_TIMERS
			timers[1][1]=clock();
			timers[1][2]+=(timers[1][1]-timers[1][0]);
			timers[3][1]=clock();
			timers[3][2]+=(timers[3][1]-timers[3][0]);
#endif // DETAILED_TIMERS
			return(-1e+200);
			//exit(0);
		      }
		    
		    F_k[k][i][j]=F_k_buff[i][j].Re();
		  }
	      
	      doubledelete(F_k_buff,num_states);
	    }
	  else
	    {
	      // Calculate F_k=Vinv*e^lambda(t(k)-t(k-1))*V :
	      for(i=0;i<num_states;i++)
		for(j=0;j<num_states;j++)
		  F_k[k][i][j]=0.0;
	  
	      // F_k=Vinv*eLambda_k*V:
	      for(i=0;i<num_states;i++)
		for(j=i%numsites2;j<num_states;j+=numsites2)
		  {
		    int start_l=num_series_feed>0 ? i%numsites2 : i;
		    for(l=start_l;l<num_states;l+=numsites2)
		      F_k[k][i][j]+=Vinv_r[i][l]*V_r[l][j]*eLambda_k_r[l];
		    // PS: F_k will also be a sparse upper triagonal matrix, just as V
		    // *if* there are no feedback loops
		  }
	    }
	}  // end if(!first)
    
      
      // Calculate x_k,(k-1) = F_k * x_(k-1),(k-1) + u_k
      for(i=0;i<num_states;i++)
	{
	  x_k_prev[i]=u_k[k][i];
	  int start_i=num_series_feed>0 ? i%numsites2 : i;
	  for(l=start_i;l<num_states;l+=numsites2)
	    x_k_prev[i]+=F_k[k][i][l]*x_k_now[(k>0 ? k-1 : 0)][l];
	}
      
      // Calculate Q_k = Vinv * Lambda_k * Vinv' :
      if(is_complex)
	{
	  for(i=0;i<num_states;i++)
	    for(j=0;j<num_states;j++)
	      {
		Qbuffer[i][j]=0.0;
		for(l=i%numsites2;l<num_states;l+=numsites2)
		  if(Vinv[i][l]!=0.0 && Lambda_k[l][j]!=0.0)
		    Qbuffer[i][j]+=Vinv[i][l]*Lambda_k[l][j];
	      }
	  for(i=0;i<num_states;i++)
	    for(j=0;j<num_states;j++)
	      {
		Q_k_buff[i][j]=0.0;
		for(l=j%numsites2;l<num_states;l+=numsites2)
		  if(Qbuffer[i][l]!=0.0 && Vinv[j][l]!=0.0)
		    Q_k_buff[i][j]+=Qbuffer[i][l]*(Vinv[j][l].conjugate());
		
		if(ABSVAL((Q_k_buff[i][j].Im()))>0.001)
		  {
#ifdef MAIN
		    static int first_cvar=1;
		    
		    cout << "Complex variance! Q_k[" << k << "][" << 
		      i << "][" << j << "]=" << 
		      Q_k_buff[i][j].Re() << "+" << 
		      Q_k_buff[i][j].Im() << "i" << endl;
		    
		    cout << endl;
		    
		    if(first_cvar)
		      {
			cout << "k=" << k << endl;
			double dt=k>0 ? t_k[k]-t_k[k-1] : 1e+200;
			cout << "dt=" << dt << endl;
			
			cout << endl;
			cout << "A:" << endl;
			for(i=0;i<num_states;i++)
			  for(j=0;j<num_states;j++)
			    cout << "A[" << i+1 << "," << j+1 << "]=" << 
			      A[i][j] << endl;
			
			cout << endl;
			cout << "lambda:" << endl;
			for(i=0;i<num_states;i++)
			  cout << "lambda[" << i+1 << "]=" << lambda[i] << endl;
			
			cout << endl;
			cout << "eLambda_k:" << endl;
			for(i=0;i<num_states;i++)
			  cout << "eLambda_k[" << i+1 << "]=" << eLambda_k[i] << endl;
			
			cout << endl;
			cout << "V:" << endl;
			for(i=0;i<num_states;i++)
			  for(j=0;j<num_states;j++)
			    cout << "V[" << i+1 << "," << j+1 << "]=" << 
			      V[i][j] << endl;
			
			cout << endl;
			cout << "Vinv:" << endl;
			for(i=0;i<num_states;i++)
			  for(j=0;j<num_states;j++)
			    cout << "Vinv[" << i+1 << "," << j+1 << "]=" << 
			      Vinv[i][j] << endl;
		    
			cout << endl;
			cout << "Sigma2:" << endl;
			for(i=0;i<num_states;i++)
			  for(j=0;j<num_states;j++)
			    cout << "Sigma2[" << i+1 << "," << j+1 << 
			      "]=" << var[i][j] << endl;
			
			cout << endl;
			cout << "Omega:" << endl;
			for(i=0;i<num_states;i++)
			  for(j=0;j<num_states;j++)
			    cout << "Omega[" << i+1 << "," << j+1 << "]=" << 
			      Omega[i][j] << endl;
			
			cout << endl;
			cout << "Lambda_k:" << endl;
			for(i=0;i<num_states;i++)
			  for(j=0;j<num_states;j++)
			    cout << "Lambda_k[" << i+1 << "," << j+1 << "]=" << 
			      Lambda_k[i][j] << endl;
			
			cout << endl;
			cout << "Qbuffer=Vinv*Lambda_k:" << endl;
			for(i=0;i<num_states;i++)
			  for(j=0;j<num_states;j++)
			    cout << "Qbuffer[" << i+1 << "," << j+1 << "]=" << 
			      Qbuffer[i][j] << endl;
			
			cout << endl;
			cout << "Q_k:" << endl;
			for(i=0;i<num_states;i++)
			  for(j=0;j<num_states;j++)
			    cout << "Q_k[" << i+1 << "," << j+1 << "]=" << 
			      Q_k_buff[i][j] << endl;
		      }    

		    first_cvar=0;
#endif // MAIN		
		    doubledelete(prior_expectancy, numobs);
		    doubledelete(resids, numobs);
		    delete [] resids_time;
#ifdef DETAILED_TIMERS
		    timers[1][1]=clock();
		    timers[1][2]+=(timers[1][1]-timers[1][0]);
		    timers[3][1]=clock();
		    timers[3][2]+=(timers[3][1]-timers[3][0]);
#endif // DETAILED_TIMERS
		    return (-1e+200);
		    
		    //exit(0);
		  }

		Q_k[k][i][j]=Q_k_buff[i][j].Re();
	      }
	}
      else
	{
	  for(i=0;i<num_states;i++)
	    for(j=0;j<num_states;j++)
	      {
		Qbuffer_r[i][j]=0.0;
		int start_l=num_series_feed>0 ? i%numsites2 : i;
		for(l=start_l;l<num_states;l+=numsites2)
		  if(Vinv_r[i][l]!=0.0 && Lambda_k_r[l][j]!=0.0)
		    Qbuffer_r[i][j]+=Vinv_r[i][l]*Lambda_k_r[l][j];
	      }
	  for(i=0;i<num_states;i++)
	    for(j=0;j<num_states;j++)
	      {
		Q_k[k][i][j]=0.0;
		int start_l=num_series_feed>0 ? j%numsites2 : j;
		for(l=start_l;l<num_states;l+=numsites2)
		  if(Qbuffer_r[i][l]!=0.0 && Vinv_r[j][l]!=0.0)
		    Q_k[k][i][j]+=Qbuffer_r[i][l]*Vinv_r[j][l];
	      }
	}
      
      // Calculate P_k_(k-1) = Q_k + F_k * P_(k-1),(k-1) * F_k' :
      for(i=0;i<num_states;i++)
	for(j=0;j<num_states;j++)
	  {
	    P_k_buffer[i][j]=0.0;
	    int start_i = num_series_feed>0 ? i%numsites2 : i;
	    for(l=start_i;l<num_states;l+=numsites2)
	      P_k_buffer[i][j]+=F_k[k][i][l]*P_k_now[(k>0 ? k-1 : 0)][l][j];
	  }
      for(i=0;i<num_states;i++)
	for(j=0;j<num_states;j++)
	  {
	    P_k_prev[k][i][j]=Q_k[k][i][j];
	    int start_j = num_series_feed>0 ? j%numsites2 : j;
	    for(l=start_j;l<num_states;l+=numsites2)
	      P_k_prev[k][i][j]+=P_k_buffer[i][l]*F_k[k][j][l];
	  }

      for(s=0;s<num_series;s++)
	if(ser[s].init_treatment &&
	   (me[k].tm==ser[s].init_time || 
	    (k==0 && me[k].tm==ser[s].meas[0].tm)))
	  {
	    int start=series_state_start[s];
	    
	    for(l=0;l<ser[s].numlayers;l++)
	      for(i=0;i<numsites;i++)
		x_k_prev[start+l*numsites+i]=ser[s].init_mu[l][i];
	    
	    for(i=0;i<ser[s].numlayers*numsites;i++)
	      for(j=0;j<ser[s].numlayers*numsites;j++)
		P_k_prev[k][start+i][start+j]=0.0;
	  }
      
      for(i=0;i<num_states;i++)
	for(j=0;j<num_states;j++)
	  if(!(P_k_prev[k][i][j]> -1e+200 && P_k_prev[k][i][j]< 1e+200))
	    {
	      doubledelete(prior_expectancy, numobs);
	      doubledelete(resids, numobs);
	      delete [] resids_time;
#ifdef DETAILED_TIMERS
	      timers[1][1]=clock();
	      timers[1][2]+=(timers[1][1]-timers[1][0]);
	      timers[3][1]=clock();
	      timers[3][2]+=(timers[3][1]-timers[3][0]);
#endif // DETAILED_TIMERS
	      return(-1e+200);
	    }

      if(debug && return_residuals)
#ifdef MAIN
	cout << "k=" << k << " me[k].num_measurements=" <<
	  me[k].num_measurements << " ret=" << ret << endl;
#else
      Rcout << "k=" << k << " me[k].num_measurements=" <<
	me[k].num_measurements << " ret=" << ret << std::endl;
#endif // MAIN
      
      if(me[k].num_measurements==1)
	{
	  unsigned int j1,j2;
	  double y_k=me[k].meanval[0]-x_k_prev[me[k].index[0]];
	  double S_k=0.0, R_k=0.0;
	  
	  // Fill out R_k (observational noise)
	  s=me[k].serie_num[0];
	  
	  if(ser[s].num_per>0)
	    for(j=0;j<ser[s].num_per;j++)
	      y_k-=ser[s].beta_sin[j]*sin(2.0*M_PI*t_k[k]/ser[s].Tper[j])-
		ser[s].beta_cos[j]*cos(2.0*M_PI*t_k[k]/ser[s].Tper[j]);
	  
	  if(me[k].sd[0]==MISSING_VALUE)
	    R_k=ser[s].obs_sd*ser[s].obs_sd;
	  else
	    R_k=me[k].sd[0]*me[k].sd[0];
	  
	  if(me[k].n[0]!=(int)MISSING_VALUE && me[k].n[0]>0)
	    R_k/=double(me[k].n[0]);
	  
	  // Calculate S_k = P_k,(k-1) + R_k
	  S_k = P_k_prev[k][me[k].index[0]][me[k].index[0]] + R_k;
	  
#ifdef MAIN
	  if(debug)
	    cout << "k=" << k << " " << me[k].meanval[0] << " " << 
	      x_k_prev[me[k].index[0]] << " " << S_k << " " << 
	      R_k << endl;

	  if(simulation_files)
	    {
	      double observation=x_k_prev[me[k].index[0]]+
		sqrt(S_k)*get_random_gauss(); // gsl_ran_ugaussian(rptr);
	      
	      if(me[k].dt!=NoHydDateTime)
		fprintf(f[s],"%04d%02d%02d/%02d%02d %d", me[k].dt.getYear(),
			me[k].dt.getMonth(),me[k].dt.getDay(),
			me[k].dt.getHour(),me[k].dt.getMinute(), me[k].site[0]);
	      else
		fprintf(f[s],"%d %9.4f ", me[k].site[0],me[k].tm);
		
	      fprintf(f[s],"%9.4f ", observation);
		
	      if(me[k].sd[0]!=MISSING_VALUE)
		fprintf(f[s], "%9.4f ", me[k].sd[0]);
	      
	      if(me[k].n[0]!=(int)MISSING_VALUE)
		fprintf(f[s],"%d ", me[k].n[0]);
	      fprintf(f[s],"\n");
	    }
#endif // MAIN

	  int s_site=me[k].site[0];
	  prior_expectancy[s*numsites+s_site][k]=x_k_prev[me[k].index[0]];
	  resids[s*numsites+s_site][k]=y_k/sqrt(S_k);
#ifndef MAIN
	  if(debug && return_residuals)
	    Rcout << "s=" << " site=" << s_site <<
	      " numsites=" << numsites << " index=" << s*numsites+s_site <<
	      " k=" << k << " y_k=" << y_k << " S_k=" << S_k <<
	      " y_k/sqrt(S_k)=" << y_k/sqrt(S_k) << std::endl;
#endif // MAIN
	      
	  // Calculate K_k = P_k,(k-1) * H_k' * inv(S_k)
	  double *K_k=new double[num_states];
	  for(j=0;j<num_states;j++)
	    K_k[j] = P_k_prev[k][j][me[k].index[0]]/S_k;
	  
	  // Calculate x_k,k = x_k,(k-1) + K_k*y_k
	  for(j=0;j<num_states;j++)
	    x_k_now[k][j] = x_k_prev[j] + K_k[j]*y_k;
	  
	  // Calculate P_k,k=P_k,(k-1) - K_k*H_k*P_k,(k-1)
	  for(j1=0;j1<num_states;j1++)
	    for(j2=0;j2<num_states;j2++)
	      P_k_now[k][j1][j2] = P_k_prev[k][j1][j2] -
		K_k[j1]*P_k_prev[k][me[k].index[0]][j2];
		
	  //bool isok=(ret>-1e+200 && ret<1e+200) ? true : false;

	  // Calculate log(f(y_k | D-1)), which is
	  // the likelihood contribution from the current measurement:
	  ret += 0.5*(log(2.0*M_PI) + log(S_k) + y_k*y_k/S_k);
	  
	  //if(!silent && isok)
	  //cout << k << " " << ret << " " << y_k << " " << S_k << endl;

	  // compensate for a log-transform in the likelihood
	  if(ser[s].pr->is_log==1) // log transformed but not pre-log transformed
	    ret -= me[k].meanval[0];

	  if(detailed_loglik)
	    {
#ifdef MAIN
	      cout << "k=" << k << " meas=" << me[k].meanval[0] <<
		" x_k_prev[" << me[k].index[0] << "]=" <<
		x_k_prev[me[k].index[0]] << " y_k=" << y_k << " S_k= " <<
		S_k << " ret= " << ret << endl;
#else
	      Rcout << "k=" << k << " meas=" << me[k].meanval[0] <<
		" x_k_prev[" << me[k].index[0] << "]=" <<
		x_k_prev[me[k].index[0]] << " y_k=" << y_k << " S_k= " <<
		S_k << " ret= " << ret << endl;
#endif
	    }
	  
#ifdef MAIN
	  if(simulation_files)
	    {
	      cout << "u_k:";
	      for(i=0;i<num_states;i++)
		cout << u_k[k][i] << " ";
	      cout << endl;

	      cout << "x_k_now:";
	      for(i=0;i<num_states;i++)
		cout << x_k_now[(k>0 ? k-1 : 0)][i] << " ";
	      cout << endl;

	      cout << "x_k_prev:";
	      for(i=0;i<num_states;i++)
		cout << x_k_prev[i] << " ";
	      cout << endl;

	      cout << k << " " << me[k].meanval[0] << " " << 
		x_k_prev[me[k].index[0]] << " " << u_k[k][me[k].index[0]] << 
		" " << ret << endl << endl;
	    }  

	  if(debug)
	    cout << "k=" << k << " l=" << -ret << endl;
#endif // MAIN	
	  
	  // Cleanup:
	  delete [] K_k;
	}
      else if(me[k].num_measurements>1)
	{
	  unsigned int n=me[k].num_measurements, i2,j1,j2;
	  double *y_k=new double[n];
	  double *zeroes_k=new double[n];
	  double **S_k=Make_matrix(n,n);
	  double **R_k=Make_matrix(n,n);

	  // Calculate y_k = Z_k - H_k*x_k,(k-1)
	  for(i=0;i<n;i++)
	    {
	      y_k[i]=me[k].meanval[i]-x_k_prev[me[k].index[i]];
	      zeroes_k[i]=0.0;
	    }

	  // Fill out R_k (observational noise)
	  for(i=0;i<n;i++)
	    {
	      s=me[k].serie_num[i];
	      
	      if(ser[s].num_per>0)
		for(j=0;j<ser[s].num_per;j++)
		  y_k[i]-=ser[s].beta_sin[j]*sin(2.0*M_PI*t_k[k]/ser[s].Tper[j])-
		    ser[s].beta_cos[j]*cos(2.0*M_PI*t_k[k]/ser[s].Tper[j]);
	      
	      if(me[k].sd[i]==MISSING_VALUE)
		R_k[i][i]=ser[s].obs_sd*ser[s].obs_sd;
	      else
		R_k[i][i]=me[k].sd[i]*me[k].sd[i];
	      
	      if(me[k].n[i]!=(int)MISSING_VALUE && me[k].n[i]>0)
		R_k[i][i]/=double(me[k].n[i]);
	    }
	  if(me[k].corr_matrix)
	    {
	      for(i=0;i<n;i++)
		for(j=0;j<n;j++)
		  if(i!=j)
		    R_k[i][j]=me[k].corr_matrix[i][j]*
		      sqrt(R_k[i][i]*R_k[j][j]);
	    }
	  
	  // Calculate S_k = P_k,(k-1) + R_k
	  for(i=0;i<n;i++)
	    for(i2=0;i2<n;i2++)
	      S_k[i][i2] = P_k_prev[k][me[k].index[i]][me[k].index[i2]] + 
		R_k[i][i2];
	  
	  double **S_k_inv=matrix_inverse(S_k, n);
	  
#ifdef MAIN
	  if(debug)
	    cout << "k=" << k << " " << me[k].meanval[0] << " " << 
	      x_k_prev[me[k].index[0]] << " " << S_k[0][0] << " " << 
	      R_k[0][0] << endl;

	  if(simulation_files)
	    {
	      double *mval=new double[n];
	      for(i=0;i<n;i++)
		mval[i]=x_k_prev[me[k].index[i]];
	      // DEBUG double **sample=sample_from_multinormal(1, mval, S_k, n, rptr);
	      double **sample=multinormal_sample(1, mval, S_k, n);
	      
	      for(i=0;i<n;i++)
		{
		  double observation=sample[0][i];

		  if(me[k].dt!=NoHydDateTime)
		    fprintf(f[s],"%04d%02d%02d/%02d%02d %d", me[k].dt.getYear(),
			    me[k].dt.getMonth(),me[k].dt.getDay(),
			    me[k].dt.getHour(),me[k].dt.getMinute(), me[k].site[0]);
		  else
		    fprintf(f[s],"%d %9.4f ", me[k].site[0],me[k].tm);
		
		  fprintf(f[s],"%9.4f ", observation);
		
		  if(me[k].sd[i]!=MISSING_VALUE)
		    fprintf(f[s], "%9.4f ", me[k].sd[i]);
		
		  if(me[k].n[i]!=(int)MISSING_VALUE)
		    fprintf(f[s],"%d ", me[k].n[i]);
		  fprintf(f[s],"\n");
		}

	      delete [] mval;
	      doubledelete(sample,1);
	    }
#endif // MAIN
	  
	  for(i=0;i<n;i++)
	    {
	      s=me[k].serie_num[i];
	      int s_site=me[k].site[i];

	      prior_expectancy[s*numsites+s_site][k]=x_k_prev[me[k].index[i]];
	      resids[s*numsites+s_site][k]=y_k[i]/sqrt(S_k[i][i]);
#ifndef MAIN
	      if(debug && return_residuals)
		Rcout << "i=" << i << " s=" << " site=" << s_site <<
		  " numsites=" << numsites << " index=" << s*numsites+s_site <<
		  " k=" << k << " y_k=" << y_k[i] << " S_k=" << S_k[i][i] <<
		  " y_k/sqrt(S_k)=" << y_k[i]/sqrt(S_k[i][i]) << std::endl;
#endif // MAIN
	    }
	  
	  // Calculate K_k = P_k,(k-1) * H_k' * inv(S_k)
	  double **K_k=Make_matrix(num_states,n);
	  for(j=0;j<num_states;j++)
	    for(i=0;i<n;i++)
	      for(i2=0;i2<n;i2++)
		K_k[j][i] += P_k_prev[k][j][me[k].index[i2]]*S_k_inv[i2][i];
	  
	  // Calculate x_k,k = x_k,(k-1) + K_k*y_k
	  for(j=0;j<num_states;j++)
	    {
	      x_k_now[k][j] = x_k_prev[j];
	      for(i=0;i<n;i++)
		x_k_now[k][j] += K_k[j][i]*y_k[i];
	    }
	  
	  // Calculate P_k,k=P_k,(k-1) - K_k*H_k*P_k,(k-1)
	  for(j1=0;j1<num_states;j1++)
	    for(j2=0;j2<num_states;j2++)
	      {
		P_k_now[k][j1][j2] = P_k_prev[k][j1][j2];
		
		for(i=0;i<n;i++)
		  P_k_now[k][j1][j2] -= K_k[j1][i]*P_k_prev[k][me[k].index[i]][j2];
	      }
	  
	  /* DEBUG stuff
	     Rcout << "loglik mainloop " << k << std::endl;

	     for(int ii=0;ii<numpar;ii++)
	     Rcout << par_name[ii] << ":" << pars[ii] << std::endl;
	  
	     show_mat_R("R_k",R_k,n,n);
	     show_mat_R("P_k_prev",P_k_prev[k],n,n);
	     show_mat_R("K_k",K_k,n,n);
	     show_mat_R("P_k_now",P_k_now[k],n,n);
	     show_mat_R("S_k",S_k,n,n);
	  */

	  if(!check_matrix(S_k,n,n))
	    {
#ifdef DETAILED_TIMERS
	      timers[1][1]=clock();
	      timers[1][2]+=(timers[1][1]-timers[1][0]);
	      timers[3][1]=clock();
	      timers[3][2]+=(timers[3][1]-timers[3][0]);
#endif // DETAILED_TIMERS
	      return -1e+200;
	    }
	  
	  // Calculate log(f(y_k | D-1)), which is
	  // the likelihood contribution from the current measurement:
	  ret -= pdf_multinormal(y_k, zeroes_k, S_k, n, true);
	  
	  for(i=0;i<n;i++)
	    {
	      s=me[k].serie_num[i];
	      // compensate for a log-transform in the likelihood
	      if(ser[s].pr->is_log==1) // log transformed but not pre-log transformed
		ret -= me[k].meanval[i];
	    }

	  if(detailed_loglik)
	    {
#ifdef MAIN
	      cout << "k=" << k << " ret= " << ret << endl;
#else
	      Rcout << "k=" << k << " ret= " << ret << endl;
#endif
	      double *x_k_prev_buffer=new double[n];
	      for(i=0;i<n;i++)
		x_k_prev_buffer[i]=x_k_prev[me[k].index[i]];
	      
	      y_k[i]=me[k].meanval[i]-x_k_prev[me[k].index[i]];
	      show_vec("meas",me[k].meanval,n);
	      show_vec("x_k_prev", x_k_prev_buffer,n);
	      show_vec("y_k", y_k,n);
	      show_mat("S_k", S_k, n, n);
	      
	      delete [] x_k_prev_buffer;
	    }
	  
#ifdef MAIN
	  if(simulation_files)
	    {
	      cout << k << " n=" << n << " m0=" << me[k].meanval[0] << 
		" ret=" << ret << endl;

	      cout << "u_k:";
	      for(i=0;i<num_states;i++)
		cout << u_k[k][i] << " ";
	      cout << endl;

	      cout << "x_k_now:";
	      for(i=0;i<num_states;i++)
		cout << x_k_now[(k>0 ? k-1 : 0)][i] << " ";
	      cout << endl;

	      cout << "F_k:";
	      for(i=0;i<num_states;i++)
		for(l=0;l<num_states;l++)
		  cout << F_k[k][i][l] << " ";
	      cout << endl;

	      cout << "indexes:";
	      for(i=0;i<n;i++)
		cout << me[k].index[i] << " ";
	      cout << endl;
	      
	      cout << "x_k_prev:";
	      for(i=0;i<n;i++)
		cout << x_k_prev[me[k].index[i]] << " ";
	      cout << endl;

	      cout << "y_k:";
	      for(i=0;i<n;i++)
		cout << y_k[i] << " ";
	      cout << endl;

	      cout << "S_k:"; 
	      for(i=0;i<n;i++)
		for(j=0;j<n;j++)
		  cout << S_k[i][j] << " ";
	      cout << endl << endl;
	    }
	  
	  if(debug)
	    cout << "k=" << k << " l=" << -ret << endl;
#endif // MAIN	  

	  // Cleanup:
	  delete [] y_k;
	  delete [] zeroes_k;
	  doubledelete(S_k,n);
	  doubledelete(S_k_inv,n);
	  doubledelete(R_k,n);
	  doubledelete(K_k, num_states);
	}
      else
	{
#ifdef MAIN
	  if(debug && t_k[k]> -10.045 && t_k[k]< -10.035)
	    for(i=0;i<num_states;i++)
	      cout << i << " " << x_k_prev[i] << " " << P_k_prev[k][i][i] << endl;
#endif // MAIN

	  //resids[k]=MISSING_VALUE;
	  for(i=0;i<num_states;i++)
	    x_k_now[k][i] = x_k_prev[i];
	  for(i=0;i<num_states;i++)
	    for(j=0;j<num_states;j++)
	      P_k_now[k][i][j] = P_k_prev[k][i][j];
	}
    }
  
#ifdef DETAILED_TIMERS
  timers[3][1]=clock();
  timers[3][2]+=(timers[3][1]-timers[3][0]);
#endif // DETAILED_TIMERS
  

  if(dosmooth)
    {
#ifdef DETAILED_TIMERS
  timers[4][0]=clock();
#endif // DETAILED_TIMERS
  

      // Smoother:
      for(i=0;i<num_states;i++)
	{
	  x_k_s[len-1][i]=x_k_now[len-1][i];
	  for(j=0;j<num_states;j++)
	    {
	      P_k_s[len-1][i][j]=P_k_now[len-1][i][j];
	      C_k[i][j]=0.0;
	    }
	}

      for(k=len-2;k>=0;k--)
	{
	  double **P_k_prevbuffer=new double*[num_states];
	  for(i=0;i<num_states;i++)
	    {
	      P_k_prevbuffer[i]=new double[num_states];
	      for(j=0;j<num_states;j++)
		P_k_prevbuffer[i][j]=P_k_prev[k+1][i][j];
	    }

	  double **P_k_plus_1_inv=matrix_inverse(P_k_prevbuffer,num_states);
	  doubledelete(P_k_prevbuffer, num_states);
	  
	  for(i=0;i<num_states;i++)
	    for(j=0;j<num_states;j++)
	      {
		PAbuffer[i][j]=0.0;
		for(l=j%numsites2;l<num_states;l+=numsites2)
		  PAbuffer[i][j]+=P_k_now[k][i][l]*F_k[k+1][j][l];
	      }
	  
	  for(i=0;i<num_states;i++)
	    for(j=0;j<num_states;j++)
	      {
		C_k[i][j]=0.0;
		for(l=0;l<num_states;l++)
		  C_k[i][j]+=PAbuffer[i][l]*P_k_plus_1_inv[l][j];
	      }
	  doubledelete(P_k_plus_1_inv,num_states);
	  
	  double *x_k_buffer=new double[num_states];
	  for(i=0;i<num_states;i++)                    
	    {
	      x_k_buffer[i]=x_k_s[k+1][i]-u_k[k+1][i];
	      for(l=0;l<num_states;l++)
		x_k_buffer[i]-=F_k[k+1][i][l]*x_k_now[k][l];
	    }
	  
	  for(i=0;i<num_states;i++)
	    {
	      x_k_s[k][i]=x_k_now[k][i];
	      for(l=0;l<num_states;l++)
		x_k_s[k][i]+=C_k[i][l]*x_k_buffer[l];
	    }
	  
	  double **P1=Make_matrix(num_states,num_states), 
	    **P2=Make_matrix(num_states,num_states);

	  for(i=0;i<num_states;i++)
	    for(j=0;j<num_states;j++)
	      {
		P1[i][j]=0.0;
		for(l=0;l<num_states;l++)
		  P1[i][j]+=P_k_s[k+1][i][l]*C_k[j][l];
	      }
	  for(i=0;i<num_states;i++)
	    for(j=0;j<num_states;j++)
	      {
		P2[i][j]=0.0;
		for(l=0;l<num_states;l++)
		  P2[i][j]+=P_k_prev[k+1][i][l]*C_k[j][l];
	      }
	  
	  for(i=0;i<num_states;i++)
	    for(j=0;j<num_states;j++)
	      {
		P_k_s[k][i][j]=P_k_now[k][i][j];

		for(l=0;l<num_states;l++)
		  P_k_s[k][i][j]+=C_k[i][l]*(P1[l][j]-P2[l][j]);
	      }

	  delete [] x_k_buffer;
	  doubledelete(P1,num_states);
	  doubledelete(P2,num_states);
	}
  
      // If periodicity in measuredments has been subtracted, add them into
      // the state again:
      for(k=0;k<len;k++)
	for(i=0;i<num_states;i++)
	  {
	    s=state_series[i];
	    
	    if(ser[s].num_per>0)
	      for(j=0;j<ser[s].num_per;j++)
		x_k_s[k][i]+=ser[s].beta_sin[j]*sin(2.0*M_PI*t_k[k]/ser[s].Tper[j])-
		  ser[s].beta_cos[j]*cos(2.0*M_PI*t_k[k]/ser[s].Tper[j]);
	  }
#ifdef DETAILED_TIMERS
  timers[4][1]=clock();
  timers[4][2]+=(timers[4][1]-timers[4][0]);
#endif // DETAILED_TIMERS
    }
  
  
  
  // sample realizations:
  
  if(do_realize)
    {
#ifdef DETAILED_TIMERS
      timers[5][0]=clock();
#endif // DETAILED_TIMERS
  
      //Rcout << "loglik realizations" << std::endl;
      
      if(!x_k_realized)
	x_k_realized=Make_matrix(meas_smooth_len, num_states);
      
      int stopped=1;
      unsigned numit=0;
      double sum_prob=0.0;
      while(stopped && numit<numit_realization)
	{
	  // sample from last place:
	  // DEBUG double **sample=sample_from_multinormal(1, x_k_s[len-1], P_k_s[len-1],
	  // num_states, rptr);
	  double **sample=multinormal_sample(1, x_k_s[len-1], P_k_s[len-1],
					     num_states);
	  
	  for(i=0;i<num_states;i++)
	    x_k_realized[len-1][i]=sample[0][i];
	  doubledelete(sample,1);
	  
	  stopped=0;

	  for(k=(len-2);k>=0 && !stopped;k--)
	    {
	      int missing=0;
	      for(i=0;i<num_states;i++)
		if(x_k_realized[k+1][i]==MISSING_VALUE)
		  missing++;
	      
	      if(missing)
		{
		  for(i=0;i<num_states;i++)
		    x_k_realized[k][i]=MISSING_VALUE;
		  stopped=1;
		}
	      else
		{
		  double sum=0.0;
		  for(i=0;i<num_states;i++)
		    for(j=0;j<num_states;j++)
		      sum+=ABSVAL((P_k_now[k][i][j]));
		  if(sum<1e-200)
		    {
		      for(i=0;i<num_states;i++)
			x_k_realized[k][i]=x_k_now[k][i];
		    }
		  else
		    {
		      double **Tk=Make_matrix(num_states, num_states), 
			*tk=new double[num_states];
		      double **Qinv=matrix_inverse(Q_k[k+1],num_states);
		      double **Pinv=matrix_inverse(P_k_now[k],num_states);
		      
		      for(i=0;i<num_states;i++)
			for(j=0;j<num_states;j++)
			  {
			    Tk[i][j]=0.0;
			    for(l=0;l<num_states;l++)
			      Tk[i][j]+=F_k[k+1][l][i]*Qinv[l][j];
			  }
		      
		      for(i=0;i<num_states;i++)
			{
			  tk[i]=0.0;
			  for(l=0;l<num_states;l++)
			    tk[i]+=Pinv[i][l]*x_k_now[k][l];
			}
		      
		      double *Tx=new double[num_states];
		      for(i=0;i<num_states;i++)
			{
			  Tx[i]=0.0;
			  for(l=0;l<num_states;l++)
			    Tx[i]+=Tk[i][l]*(x_k_realized[k+1][l]-u_k[k+1][l]);
			}
		      
		      if(t_k[k]>=t_k[k+1])
			{
#ifdef MAIN
			  cerr << k << " " << Tx[2] << " " <<
			    x_k_realized[k+1][2] << " t_k[k]=" << t_k[k] <<
			    " t_k[k+1]" << t_k[k+1] << " t_k[k+1]-t_k[k]=" <<
			    t_k[k+1]-t_k[k] << endl;
			  exit(0);
#else
			  Rcout << k << " " << Tx[2] << " " <<
			    x_k_realized[k+1][2] << " t_k[k]=" << t_k[k] <<
			    " t_k[k+1]" << t_k[k+1] << " t_k[k+1]-t_k[k]=" <<
			    t_k[k+1]-t_k[k] << std::endl;
#endif // MAIN
			}
		      
		      double **Sk_inv=new double*[num_states];
		      for(i=0;i<num_states;i++)
			{
			  Sk_inv[i]=new double[num_states];
			  for(j=0;j<num_states;j++)
			    {
			      Sk_inv[i][j]=Pinv[i][j];
			      for(l=0;l<num_states;l++)
				Sk_inv[i][j]+=Tk[i][l]*F_k[k+1][l][j];
			    }
			}
		      double **Sk=matrix_inverse(Sk_inv,num_states);
		      
		      double *raw_mean=new double[num_states],
			*mean=new double[num_states];
		      for(i=0;i<num_states;i++)
			raw_mean[i]=tk[i]+Tx[i];
		      
		      for(i=0;i<num_states;i++)
			{
			  mean[i]=0.0;
			  for(l=0;l<num_states;l++)
			    mean[i]+=Sk[i][l]*raw_mean[l];
			}

		      if(!check_matrix(Sk,num_states,num_states))
			{
#ifdef DETAILED_TIMERS
			  timers[1][1]=clock();
			  timers[1][2]+=(timers[1][1]-timers[1][0]);
			  timers[5][1]=clock();
			  timers[5][2]+=(timers[5][1]-timers[5][0]);
#endif // DETAILED_TIMERS
			  return -1e+200;
			}
		      
		      double *e=double_eigenvalues(Sk, num_states);
		      
		      if(!e)
			{
#ifdef MAIN
			  cout << t_k[k] << " - " << t_k[k+1] << endl;
			  cout << "Couldn't evaluate eigenvalues!" << endl;
			  for(i=0;i<num_states;i++)
			    {
			      for(j=0;j<num_states;j++)
				cout << P_k_now[k][i][j] << " ";
			      cout << endl;
			    }
#endif // MAIN
			  for(i=0;i<num_states;i++)
			    x_k_realized[k][i]=MISSING_VALUE;
			}
		      else 
			{
			  int neg_eigen=0;
			  for(i=0;i<num_states && !stopped;i++)
			    if(e[i]<=0.0)
			      {
#ifdef MAIN
				cout << t_k[k] << " - " << t_k[k+1] << endl;
				cout << "Eigenvalue " << i+1 << " = " << 
				  e[i] << endl;
				cout << "k=" << k << endl;
#endif // MAIN				

				for(i=0;i<num_states;i++)
				  x_k_realized[k][i]=MISSING_VALUE;
				neg_eigen=1;
			      }
			  
			  if(e)
			    delete [] e;
			  
			  if(!neg_eigen)
			    {
			      if(!check_matrix(Sk,num_states,num_states))
				{
#ifdef DETAILED_TIMERS
				  timers[1][1]=clock();
				  timers[1][2]+=(timers[1][1]-timers[1][0]);
				  timers[5][1]=clock();
				  timers[5][2]+=(timers[5][1]-timers[5][0]);
#endif // DETAILED_TIMERS
				  return -1e+200;
				}
			      
			      double det=matrix_determinant(Sk, num_states);
			      if(det<1e-200)
				{
#ifdef MAIN
				  cout << t_k[k] << " - " << t_k[k+1] << endl;
				  cout << "det=" << det << endl;
#endif // MAIN
				  for(i=0;i<num_states;i++)
				    x_k_realized[k][i]=MISSING_VALUE;
				}
			      else
				{
				  //sample=
				  //sample_from_multinormal(1, mean, Sk, 
				  //num_states, rptr);
				  sample=
				    multinormal_sample(1, mean, Sk, 
						       num_states);
				  for(i=0;i<num_states;i++)
				    x_k_realized[k][i]=sample[0][i];
				  doubledelete(sample,1);
				}
			    }
			}

		      doubledelete(Qinv,num_states);
		      doubledelete(Pinv,num_states);
		      doubledelete(Sk,num_states);
		      doubledelete(Tk,num_states);
		      doubledelete(Sk_inv,num_states);
		      delete [] Tx;
		      delete [] tk;
		      delete [] raw_mean;
		      delete [] mean;
		    }
		}
	    }

	  if(real_strat!=NO_CENSORING && !stopped)
	    {
	      unsigned int *k2=new unsigned int[numsites];
	      // stores the index of the last site specific measurement
	      for(i=0;i<numsites;i++)
		k2[i]=len-1;
	      
	      for(k=(len-2);k>=0 && !stopped;k--)
		{
		  // Find the index of the last site specific measurement
		  // before the currently examined time point:
		  for(i=0;i<numsites && !stopped;i++)
		    {
		      while(k2[i]>=0 && 
			    (k2[i]>=k || 
			     (me[k2[i]].site!=NULL &&
			      me[k2[i]].site[0]!=(int)i) || 
			     (me[k2[i]].meanval!=NULL && 
			      me[k2[i]].meanval[0]==MISSING_VALUE)))
			k2[i]--;
		      
		      if(k2[i]>=0)
			{
			  if(real_strat==ASCENDING_CENSORING &&
			     x_k_realized[k][i]>=x_k_realized[k2[i]][i])
			    {
			      stopped=1;
#ifdef MAIN
			      if(!silent)
				printf("iteration %d/%d - site %d: "
				       "(t=%8.3f,x=%8.3f)>=(t=%8.3f,x=%8.3f)\n",
				       numit, numit_realization,i, 
				       t_k[k], x_k_realized[k][i],
				       t_k[k2[i]],x_k_realized[k2[i]][i]);
#endif  // MAIN
			    }
				
			  if(real_strat==DESCENDING_CENSORING &&
			     x_k_realized[k][i]<=x_k_realized[k2[i]][i])
			    {
			      stopped=1;
#ifdef MAIN
			      if(!silent)
				printf("iteration %d/%d - site %d: "
				       "(t=%8.3f,x=%8.3f)<=(t=%8.3f,x=%8.3f)\n",
				       numit, numit_realization,i, 
				       t_k[k], x_k_realized[k][i],
				       t_k[k2[i]],x_k_realized[k2[i]][i]);
#endif  // MAIN
			    }
			}
		    }
		}

	      sum_prob+=1.0/double((numit+1)*(numit+2));

	      delete [] k2;
	    }

	  if(stopped)
	    numit++;
	}

      if(stopped || numit>=numit_realization)
	{
	  doubledelete(x_k_realized,num_states);
	  x_k_realized=NULL;
	  ret -= log(1.0-sum_prob);
	}
      else
	ret += log(double(numit+1))+log(double(numit+2));
      
#ifdef DETAILED_TIMERS
      timers[5][1]=clock();
      timers[5][2]+=(timers[5][1]-timers[5][0]);
#endif // DETAILED_TIMERS
    }

  if(residual_analysis)
    analyze_residuals(resids[0], t_k, dt_k,len, res_filestart);
  if(return_residuals && residuals!=NULL &&
     resid_len!=NULL && residuals_time!=NULL)
    {
      *residuals_time=resids_time;
      *prior_expected_values=prior_expectancy;
      *residuals=resids;
      *resid_numcolumns=numobs;
      *resid_len=len;
    }
  else
    {
      doubledelete(prior_expectancy, numobs);
      doubledelete(resids, numobs);
      delete [] resids_time;
    }
  
  if(dt_k)
    delete [] dt_k;

  if(dosmooth)
    {
      keep_x_and_P(t_k, x_k_s, P_k_s,len);
    }

  // Cleanup:
  if(is_complex)
    {
      if(lambda)
	delete [] lambda;
      doubledelete(V,num_states);
      doubledelete(Vinv,num_states);
      doubledelete(VinvLambdaV_A,num_states);
      doubledelete(Omega,num_states);
      doubledelete(Vvar,num_states);
      doubledelete(Qbuffer,num_states);
      doubledelete(Q_k_buff,num_states);
      doubledelete(Lambda_k,num_states);
      if(eLambda_k)
	delete [] eLambda_k;
      if(W)
	delete [] W;
      if(u_k_buff)
	delete [] u_k_buff;
    }
  else
    {
      delete [] lambda_r;
      doubledelete(V_r,num_states);
      doubledelete(Vinv_r,num_states);
      doubledelete(Omega_r,num_states);
      doubledelete(Vvar_r,num_states);
      doubledelete(Qbuffer_r,num_states);
      doubledelete(Lambda_k_r,num_states);
      if(eLambda_k_r)
	delete [] eLambda_k_r;
      if(W_r)
	delete [] W_r;
      if(u_k_buff_r)
	delete [] u_k_buff_r;
    }
  
  delete [] t_k;
  doubledelete(var,num_states);
  doubledelete(A,num_states);
  doubledelete(P_k_buffer,num_states);
  
  for(k=0;k<len;k++)
    {
      doubledelete(P_k_now[k],num_states);
      doubledelete(P_k_prev[k],num_states);
      doubledelete(F_k[k],num_states);
      doubledelete(Q_k[k],num_states);
      if(dosmooth)
	doubledelete(P_k_s[k],num_states);
    }
  delete [] Q_k;
  delete [] P_k_now;
  delete [] P_k_prev;
  delete [] F_k;
  doubledelete(x_k_now,len);
  doubledelete(u_k,len);
  if(dosmooth)
    {
      doubledelete(x_k_s,len);
      delete [] P_k_s;
      doubledelete(C_k,num_states);
      doubledelete(PAbuffer,num_states);
    }
  
  delete [] x_k_prev;
  delete [] u_k;
  delete [] m;

#ifdef MAIN
  if(f)
    {
      for(s=0;s<num_series;s++)
	fclose(f[s]);
      delete [] f;
      f=NULL;
    }
#endif // MAIN
  
  if(talkative_likelihood)
#ifdef MAIN
    cout << "ll=" << -ret << endl;
#else  
  Rcout << "ll=" << -ret << std::endl;
#endif // MAIN
  
#ifdef DETAILED_TIMERS
  timers[1][1]=clock();
  timers[1][2]+=(timers[1][1]-timers[1][0]);
#endif // DETAILED_TIMERS
  
  // return the loglikelihood:
  return -ret;
}



double minusloglik(double *pars)
{
  return -loglik(pars);
}



// **************************************************************
// logprob: Returns log-probability (prior times likelihood).
// Input: transformed parameter vector,  
// and 'temperature' (used for tempering purposes).
// **************************************************************

double logprob(params &par,double T, int dosmooth, int do_realize, int doprint)
{
  unsigned int s,i,l;
  
  // fetch the log-likelihood:
  double ll=loglik(par.param, dosmooth, do_realize), prp=0.0;;

  // sanity check on the likelihood:
  if(!(ll>-1e+199 && ll<1e+199))
    return -1e+200;


  
#ifdef DETAILED_TIMERS
  timers[9][0]=clock();
#endif // DETAILED_TIMERS
  

  // handle the pull priors separately, according to the identification
  // prior handling strategy:

  for(s=0;s<num_series;s++)
    {
      prior *pr=ser[s].pr;
      for(l=0;l<ser[s].numlayers;l++)
	if(!ser[s].time_integral[l] && 
	   !(l==(ser[s].numlayers-1) && ser[s].no_pull_lower))
	  {
	    switch(id_strategy)
	      {
	      case ID_NONE:
		{
		  if(!ser[s].allow_positive_pulls)
		    {
		      if(l==(ser[s].numlayers-1) && ser[s].no_pull_lower) 
			{
			  // NOP
			}
		      else if(ser[s].regional_pull[l])
			{
			  for(i=0;i<numsites;i++)
			    {
			      double ldt=-log(ser[s].pull[l][i]);
			      prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
				0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/
				pr->ldt_s/pr->ldt_s;
			    }
			}
		      else if(ser[s].indicator_pull[l])
			{
			  int done0=0,done1=0;
			  for(i=0;i<numsites;i++)
			    {
			      if((indicator_array[i]==0 && !done0) ||
				 (indicator_array[i]==1 && !done1))
				{
				  double ldt=-log(ser[s].pull[l][i]);
				  prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
				    0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/
				    pr->ldt_s/pr->ldt_s;
				  
				  if(indicator_array[i]==0)
				    done0=1;
				  else if(indicator_array[i]==1)
				    done1=1;
				}
			    }
			}
		      else
			{
			  double ldt=-log(ser[s].pull[l][0]);
			  prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
			    0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/
			    pr->ldt_s/pr->ldt_s;
			}
		    }
		  else
		    {
		      double pull_mu=0.0, pull_sd=1.0/pr->dt_1/1.96;
		      
		      if(l==(ser[s].numlayers-1) && ser[s].no_pull_lower) 
			{
			  // NOP
			}
		      else if(ser[s].regional_pull[l])
			{
			  for(i=0;i<numsites;i++)
			    {
			      double pull=ser[s].pull[l][i];
			      prp += -0.5*log(2.0*M_PI) - log(pull_sd) -
				0.5*(pull - pull_mu)*(pull - pull_mu)/
				pull_sd/pull_sd;
			    }
			}
		      else if(ser[s].indicator_pull[l])
			{
			  int done0=0,done1=0;
			  for(i=0;i<numsites;i++)
			    {
			      if((indicator_array[i]==0 && !done0) ||
				 (indicator_array[i]==1 && !done1))
				{
				  double pull=ser[s].pull[l][i];
				  prp += -0.5*log(2.0*M_PI) - log(pull_sd) -
				    0.5*(pull - pull_mu)*(pull - pull_mu)/
				    pull_sd/pull_sd;

				  if(indicator_array[i]==0)
				    done0=1;
				  else if(indicator_array[i]==1)
				    done1=1;
				}
			    }
			}
		      else
			{
			  double pull=ser[s].pull[l][0];
			  prp += -0.5*log(2.0*M_PI) - log(pull_sd) -
			    0.5*(pull - pull_mu)*(pull - pull_mu)/
			    pull_sd/pull_sd;
			}
		    }
		  break;
		}
	      case ID_CUT_LOWER:
	      case ID_SUB_LOWER:
		{
		  if(l==(ser[s].numlayers-1) && ser[s].no_pull_lower) 
		    ; // NOP
		  else if(ser[s].regional_pull[l])
		    {
		      for(i=0;i<numsites;i++)
			{
			  double ldt=-log(ser[s].pull[l][i]);
			  
			  if(l==(ser[s].numlayers-1))
			    prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
			      0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/
			      pr->ldt_s/pr->ldt_s;
			  else
			    {
			      double ldt_old=-log(ser[s].pull[l+1][i]);
			      
			      if(ldt > ldt_old)
				{
#ifdef DETAILED_TIMERS
				  timers[9][1]=clock();
				  timers[9][2]+=(timers[9][1]-timers[9][0]);
#endif // DETAILED_TIMERS
				  return -1e+200;
				}
			      
			      if(id_strategy==ID_CUT_LOWER)
				prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
				  0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/
				  pr->ldt_s/pr->ldt_s -
				  log(standard_normal_cdf((ldt_old-pr->ldt_m)/
							  pr->ldt_s));
			      //log(gsl_cdf_ugaussian_P((ldt_old-pr->ldt_m)/
			      //pr->ldt_s));
			      else if(id_strategy==ID_SUB_LOWER)
				prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
				  log(exp(ldt_old-ldt)-1.0)-
				  0.5*(log(exp(ldt_old)-exp(ldt)) - 
				       pr->ldt_m)*
				  (log(exp(ldt_old)-exp(ldt)) - pr->ldt_m)/
				  pr->ldt_s/pr->ldt_s;
			    }
			}
		    }
		  else if(ser[s].indicator_pull[l])
		    {
		      int done0=0,done1=0;
		      for(i=0;i<numsites;i++)
			{
			  if((indicator_array[i]==0 && !done0) ||
			     (indicator_array[i]==1 && !done1))
			    {
			      double ldt=-log(ser[s].pull[l][i]);
			      
			      if(l==(ser[s].numlayers-1))
				prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
				  0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/
				  pr->ldt_s/pr->ldt_s;
			      else
				{
				  double ldt_old=-log(ser[s].pull[l+1][i]);
				  
				  for(unsigned int i2=0;i2<numsites;i2++)
				    if(i2!=i && indicator_array[i]==
				       indicator_array[i2])
				      ldt_old=MINIM(ldt_old,-log(ser[s].pull[l+1][i]));
				  
				  if(ldt > ldt_old)
				    {
#ifdef DETAILED_TIMERS
				      timers[9][1]=clock();
				      timers[9][2]+=(timers[9][1]-timers[9][0]);
#endif // DETAILED_TIMERS    
				      return -1e+200;
				    }
				  
				  if(id_strategy==ID_CUT_LOWER)
				    prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
				      0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/
				      pr->ldt_s/pr->ldt_s -
				      log(standard_normal_cdf((ldt_old - pr->ldt_m)/pr->ldt_s));
				  // log(gsl_cdf_ugaussian_P((ldt_old - pr->ldt_m)/pr->ldt_s));
				  else if(id_strategy==ID_SUB_LOWER)
				    prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
				      log(exp(ldt_old-ldt)-1.0)-
				      0.5*(log(exp(ldt_old)-exp(ldt)) - 
					   pr->ldt_m)*
				      (log(exp(ldt_old)-exp(ldt)) - pr->ldt_m)/
				      pr->ldt_s/pr->ldt_s;
				}
				
			      if(indicator_array[i]==0)
				done0=1;
			      else if(indicator_array[i]==1)
				done1=1;
			    }
			}
		    }
		  else
		    {
		      double ldt=-log(ser[s].pull[l][0]);
		      
		      if(l==(ser[s].numlayers-1))
			prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
			  0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/pr->ldt_s/pr->ldt_s;
		      else 
			{
                          bool allzero_nextlayer=true;
			  for(i=0;i<numsites;i++)
			    if(ser[s].pull[l+1][i]>0.0)
			      allzero_nextlayer=false;
			  
			  if(!allzero_nextlayer)
			    {
			      double ldt_old=-log(ser[s].pull[l+1][0]+1e-100);
			      
			      for(i=1;i<numsites;i++)
				ldt_old=MINIM(ldt_old,-log(ser[s].pull[l+1][i]+1e-100));
			      
			      if(ldt > ldt_old)
				{
#ifdef DETAILED_TIMERS
				  timers[9][1]=clock();
				  timers[9][2]+=(timers[9][1]-timers[9][0]);
#endif // DETAILED_TIMERS
				  return -1e+200;
				}
			      
			      if(id_strategy==ID_CUT_LOWER)
				prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
				  0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/
				  pr->ldt_s/pr->ldt_s -
				  log(standard_normal_cdf((ldt_old - pr->ldt_m)/
							  pr->ldt_s));
			      //log(gsl_cdf_ugaussian_P((ldt_old - pr->ldt_m)/
			      //pr->ldt_s));
			      else if(id_strategy==ID_SUB_LOWER)
				prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
				  log(exp(ldt_old-ldt)-1.0)-
				  0.5*(log(exp(ldt_old)-exp(ldt)) - pr->ldt_m)*
				  (log(exp(ldt_old)-exp(ldt)) - pr->ldt_m)/
				  pr->ldt_s/pr->ldt_s;
			    }
			  else
			    {
			      prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
				0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/
				pr->ldt_s/pr->ldt_s;
			    }
			}
		    }
		  
		  break;
		}
	      case ID_CUT_UPPER:
	      case ID_ADD_UPPER:
		{
		  if(l==(ser[s].numlayers-1) && ser[s].no_pull_lower) 
		    ; // NOP
		  else if(ser[s].regional_pull[l])
		    {
		      for(i=0;i<numsites;i++)
			{
			  double ldt=-log(ser[s].pull[l][i]);
			  
			  if(l==0)
			    prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
			      0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/
			      pr->ldt_s/pr->ldt_s;
			  else
			    {
			      double ldt_old=-log(ser[s].pull[l-1][i]);
			      
			      if(ldt < ldt_old)
				{
#ifdef DETAILED_TIMERS
				  timers[9][1]=clock();
				  timers[9][2]+=(timers[9][1]-timers[9][0]);
#endif // DETAILED_TIMERS
				  return -1e+200;
				}
			      
			      if(id_strategy==ID_CUT_UPPER)
				prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
				  0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/
				  pr->ldt_s/pr->ldt_s -
				  log(1.0-standard_normal_cdf((ldt_old - pr->ldt_m)/pr->ldt_s));
			      // log(gsl_cdf_ugaussian_Q((ldt_old - pr->ldt_m)/pr->ldt_s));
			      else if(id_strategy==ID_ADD_UPPER)
				prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
				  log(1.0-exp(ldt_old-ldt))-
				  0.5*(log(exp(ldt)-exp(ldt_old)) - pr->ldt_m)*
				  (log(exp(ldt)-exp(ldt_old)) - pr->ldt_m)/
				  pr->ldt_s/pr->ldt_s;
			    }
			}
		    }
		  else if(ser[s].indicator_pull[l])
		    {
		      int done0=0,done1=0;
		      for(i=0;i<numsites;i++)
			{
			  if((indicator_array[i]==0 && !done0) ||
			     (indicator_array[i]==1 && !done1))
			    {
			      double ldt=-log(ser[s].pull[l][i]);
			      
			      if(l==0)
				prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
				  0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/
				  pr->ldt_s/pr->ldt_s;
			      else
				{
				  double ldt_old=-log(ser[s].pull[l-1][i]);
				  
				  for(unsigned int i2=0;i2<numsites;i2++)
				    if(i2!=i && 
				       indicator_array[i]==indicator_array[i2])
				      ldt_old=MAXIM(ldt_old,
						    -log(ser[s].pull[l-1][i]));
				  
				  if(ldt < ldt_old)
				    {
#ifdef DETAILED_TIMERS
				      timers[9][1]=clock();
				      timers[9][2]+=(timers[9][1]-timers[9][0]);
#endif // DETAILED_TIMERS
				      return -1e+200;
				    }
				  
				  if(id_strategy==ID_CUT_UPPER)
				    prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
				      0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/
				      pr->ldt_s/pr->ldt_s -
				      log(1.0-standard_normal_cdf((ldt_old - pr->ldt_m)/pr->ldt_s));
				  // log(gsl_cdf_ugaussian_Q((ldt_old - pr->ldt_m)/pr->ldt_s));
				  else if(id_strategy==ID_ADD_UPPER)
				    prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
				      log(1.0-exp(ldt_old-ldt))-
				      0.5*(log(exp(ldt)-exp(ldt_old)) - 
					   pr->ldt_m)*
				      (log(exp(ldt)-exp(ldt_old)) - pr->ldt_m)/
				      pr->ldt_s/pr->ldt_s;
				}
			      
			      if(indicator_array[i]==0)
				done0=1;
			      else if(indicator_array[i]==1)
				done1=1;
			    }
			}
		    }
		  else
		    {
		      double ldt=-log(ser[s].pull[l][0]);
		      
		      if(l==0)
			prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
			  0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/
			  pr->ldt_s/pr->ldt_s;
		      else
			{
			  double ldt_old=-log(ser[s].pull[l-1][0]);
			  
			  for(i=1;i<numsites;i++)
			    ldt_old=MAXIM(ldt_old,-log(ser[s].pull[l-1][i]));
			  
			  if(ldt < ldt_old)
			    {
#ifdef DETAILED_TIMERS
			      timers[9][1]=clock();
			      timers[9][2]+=(timers[9][1]-timers[9][0]);
#endif // DETAILED_TIMERS
			      return -1e+200;
			    }
			  
			  if(id_strategy==ID_CUT_UPPER)
			    prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
			      0.5*(ldt - pr->ldt_m)*(ldt - pr->ldt_m)/
			      pr->ldt_s/pr->ldt_s -
			      log(1.0-standard_normal_cdf((ldt_old - pr->ldt_m)/pr->ldt_s));
			  // log(gsl_cdf_ugaussian_Q((ldt_old - pr->ldt_m)/pr->ldt_s));
			  else if(id_strategy==ID_ADD_UPPER)
			    prp += -0.5*log(2.0*M_PI) - log(pr->ldt_s) -
			      log(1.0-exp(ldt_old-ldt))-
			      0.5*(log(exp(ldt)-exp(ldt_old)) - pr->ldt_m)*
			      (log(exp(ldt)-exp(ldt_old)) - pr->ldt_m)/
			      pr->ldt_s/pr->ldt_s;
			}
		    }
		  
		  break;
		}
		
	      default:
#ifdef MAIN
		cerr << "Unknown identification-prior strategy!" << endl;
		exit(0);
#else
		Rcout << "Unknown identification-prior strategy!" << std::endl;
#endif // MAIN
		break;
	      }
	  }
    }

  if(some_pairwise_correlations)
    {
      for(s=0;s<num_series;s++)
	for(l=0;l<ser[s].numlayers;l++)
	  if(ser[s].sigma_pairwise_correlated[l])
	    prp-=log(p_pos_site_sigma2);
    }

  if(num_series_corr>0)
    prp-=log(p_pos_series_sigma2);
  
  // prior for all non-pull parameters::
  for(i=0;i<numpar && prp>-1e+200;i++) // traverse the parameters:
    {
      s=par_series[i];
      prior *pr=ser[s].pr;
      
      // Fetch the current transformed parameter value (which should be
      // normally distributed in the prior):
      double curr_val=par.param[i]; 
      
      // depending on the parameter type, update the prior density:
      switch(par_type[i])
	{
	case MU:
	  {
	    prp += -0.5*log(2.0*M_PI) - log(pr->mu_s) -
	      0.5*(curr_val-pr->mu_m)*(curr_val-pr->mu_m)/pr->mu_s/pr->mu_s;
	    break;
	  }
	case TRIG: // Use prior uncertainty as standard deviation plus
		   // expectancy of MU
	  {
	    prp += -0.5*log(2.0*M_PI) - log(pr->mu_s+pr->mu_m) -
	      0.5*curr_val*curr_val/(pr->mu_s+pr->mu_m)/(pr->mu_s+pr->mu_m);
	    break;
	  }
	case LIN_T:
	  {
	    prp += -0.5*log(2.0*M_PI) - log(pr->lin_s) -
	      0.5*(curr_val-pr->lin_m)*(curr_val-pr->lin_m)/pr->lin_s/pr->lin_s;
	    break;
	  }
	case DT:
	  {
	    break;
	  }
	case SIGMA:
	  {
	    prp+=  -0.5*log(2.0*M_PI) - log(pr->ls_s) -
	      0.5*(curr_val-pr->ls_m)*(curr_val-pr->ls_m)/pr->ls_s/pr->ls_s;
	    break;
	  }
	case OBS_SD:
	  {
	    prp+=  -0.5*log(2.0*M_PI) - log(pr->los_s) -
	      0.5*(curr_val - pr->los_m)*(curr_val - pr->los_m)/pr->ls_s/pr->los_s;
	    //cout << curr_val << " " << pr->los_m << " " << pr->los_s << endl;
	    break;
	  }
	case CORR:
	  {
	    prp += -0.5*log(2.0*M_PI) - log(2.0) -
	      0.5*curr_val*curr_val/2.0/2.0;
	    break;
	  }
	case PAIR_CORR:
	  {
	    prp += -0.5*log(2.0*M_PI) - log(2.0) -
	      0.5*curr_val*curr_val/2.0/2.0;
	    break;
	  }
	case SERIES_CORR:
	  {
	    prp += -0.5*log(2.0*M_PI) - log(2.0) -
	      0.5*curr_val*curr_val/2.0/2.0;
	    break;
	  }
	case BETA:
	  {
	    prp += -0.5*log(2.0*M_PI) - log(pr->beta_s) -
	      0.5*(curr_val-pr->beta_m)*(curr_val-pr->beta_m)/pr->beta_s/pr->beta_s;
	    break;
	  }
	case INDICATOR:
	  {
	    prp += log(0.5);
	    break;
	  }
	case INIT:
	  {
	    prp +=  -0.5*log(2.0*M_PI) - log(pr->init_s) -
	      0.5*(curr_val-pr->init_m)*(curr_val-pr->init_m)/pr->init_s/pr->init_s;
	    break;
	  }
	}
    }

  
#ifdef DETAILED_TIMERS
  timers[9][1]=clock();
  timers[9][2]+=(timers[9][1]-timers[9][0]);
#endif // DETAILED_TIMERS
  

  // sanity check on the prior:
  if(prp<=-1e+200)
    return -1e+200;
  // sanity check on the likelihood:
  if(ll<=-1e+200)
    return -1e+200;
  
  // add logarithmic prior and likelihood and divide by
  // the tempering temperature:
  double ret=(ll+prp)/T;

  if(doprint)
#ifdef MAIN
    cout << "ll=" << ll << " prp=" << prp << endl;
#else
  Rcout << "ll=" << ll << " prp=" << prp << std::endl;
#endif // MAIN

  // update the parameter likelihood and probability density:
  par.log_lik=ll;
  par.log_prob=ll+prp;

  // return the tempered log(prior*likelihood):
  return ret;
}



// ******************************************************
// newsample: Sample one new parameter set, using 
// Random Walk Metropolis and parallell tempering on
// the previous sample.
// sample: starts as the previous sample, ends as the
//         current sample. It's an array, since you can
//         have several tempering chains.
// logprob: start as previous logprobs 
//          (log(prior*likelihood) for each chain)
//          and ends as the current logprobs.
// acc: acceptance indicator array. Is set to one for 
//      each parameter that has it's new value accepted.
//      Used for the adaptive phase of the burn-in.
// rw: Random Walk standard deviations for the proposal
//     density of the Metropolis algorithm.
// T: temperature array
// numtemp: Number of tempering chains
// swaps: Incremented each time a swapping of tempering
//        chains takes place.
// ******************************************************
void newsample(params *sample, double *log_prob, params *acc, params *rw,
	       double *T, unsigned int numtemp, int *swaps,
	       int dosmooth, int do_realize)
{
  if(numtemp>1 && drand()<0.1) // tempering swap?
    {
      // choose a random chain below the maximum:
      int index=(int) floor(double(numtemp-1)*drand());
      // Find the logprob when a switch takes place:
      double p1=logprob(sample[index], T[index+1]);
      double p2=logprob(sample[index+1], T[index]);
      
      // Acceptance-rejection term:
      if(log(drand())<p1+p2-log_prob[index]-log_prob[index+1])
	{
	  // Swap chains:
	  log_prob[index+1]=p1;
	  log_prob[index]=p2;
	  params buffer(sample+index);
	  sample[index].copy(sample+index+1);
	  sample[index+1].copy(&buffer);
	  swaps[index]++;
	}
    }
  else // normal RW Metropolis sampling:
    for(unsigned int t=0;t<numtemp;t++) // traverse the tempering chains
      {
	// Make a new parameter structure from the old:
	params newparams(sample+t);
	double new_logprob;
	double sT=sqrt(T[t]); // tempering contribution to the 
	// random walk standard deviation
	
	for(unsigned int i=0;i<numpar;i++) // traverse the parameters
	  {
	    // if binary, try to switch the indicator, else
	    // sample from the normal distribution with the
	    // stated standard deviation found in 'rw':
	    if(par_trans_type[i]==T_BINARY)
	      {
		if(sample[t].param[i]!=0.0)
		  newparams.param[i]=0.0;
		else
		  newparams.param[i]=1.0;
	      }
	    else
	      newparams.param[i]=sample[t].param[i]+sT*rw->param[i]*get_random_gauss();

	    // Update the logprob (log(prior*likelihood)):
	    new_logprob=logprob(newparams, T[t]);

	    // Acceptance/rejection term:
	    if(log(drand())<new_logprob-log_prob[t])
	      {
		// copy the new value over to the 'sample' structure:
		sample[t].copy(&newparams);
		// update the logprob:
		log_prob[t]=new_logprob;
		// if this is the lowest chain (the chain of the
		// model we are interested in), then set the
		// acceptances indicator for this parameter:
		if(t==0)
		  acc->param[i]=1.0;
	      }
	    else // rejection
	      {
		// if this is the lowest chain (the chain of the
		// model we are interested in), then nullify the
		// acceptances indicator for this parameter:
		if(t==0)
		  acc->param[i]=0.0;
		// put the old sample value into the 
		// new sample structure
		newparams.param[i]=sample[t].param[i];
	      }
	  }

	// update the logprob:
	sample[t].log_prob=log_prob[t];
      }

  if(dosmooth)
    logprob(sample[0], T[0], 1, do_realize);
  else
    logprob(sample[0], T[0]);
}

double init_par(int i) // i=parameter number
{
  double ret=0.0;
  
  switch(par_type[i])
    {
      // fetch parameter suggestions according to 
      // the parameter type:
    case INDICATOR:
      ret=floor(2.0*drand());
      break;
    case MU:
    case LIN_T:
    case TRIG:
      ret=-3.0+6.0*drand();
      break;
    case DT:
    case SIGMA:
    case OBS_SD:
    case CORR:
    case INIT:
      ret=-4.0+8.0*drand();
      break;
    case PAIR_CORR:
    case SERIES_CORR:
      ret=-0.1+0.2*drand();
      break;
    case BETA:
      ret=-1.0+2.0*drand();
    }

  return ret;
}

// ***************************************************
// mcmc: Perform Markov Chain Monte Carlo sampling.
// numsamples: number of samples.
// burnin: burn-in period.
// indep: number of MCMC iterations between each
//        sample fetched.
// numtemp: Number of tempering chains.
// ****************************************************
params *layer_mcmc(unsigned int numsamples, unsigned int burnin,
		   unsigned int indep,unsigned int numtemp,
		   bool do_importance,
		   int dosmooth, int do_realization,  
		   char *realization_file_start, double ****x_, 
		   double *startpar, double T_ground,
		   double *model_loglik, 
		   double *model_dic1,double *eff_num_param1, 
		   double *model_dic2,double *eff_num_param2)
{
  loglik(NULL); // gives us the global 'numpar', 
  // number of parameters
  
  // Return array:
  params *ret=new params[numsamples];
  
  double *ll=new double[numsamples];
  double *lp=new double[numsamples];
  
  // Current set of samples (one for each tempering chain) and
  // acceptance indicators and random walk standard deviations:
  params *pars=new params[numtemp], acc(numpar), 
    totalacc(numpar), rw(numpar);
  
  // logprob, i.e. log(prior*likelihood):
  double *log_prob=new double[numtemp];
  
  // indexes:
  unsigned int i,j,k;
  unsigned int t;
  
  // Temperatures for each tempering chain:
  double *T=new double[numtemp];
  
  // Number of swaps between one chain and the next:
  int *swaps=new int[numtemp];
  
  // Initialize GSL random generator:
  //gsl_rng *rptr=gsl_rng_alloc(gsl_rng_rand48);
  //gsl_rng_set(rptr, rand()); 
  
  // traverse the tempering chains:
  for(t=0;t<numtemp;t++)
    {
      // T(t)=T_ground^(t-1), T_ground is default 2
      T[t]=pow(T_ground,double(t));
      swaps[t]=0;
    }
  
  // initialize the random walk standard deviations:
  for(i=0;i<numpar;i++)
    rw.param[i]=1.0;
  
  
#ifdef DETAILED_TIMERS
  timers[20][0]=clock();
#endif // DETAILED_TIMERS
  
  
  // **************************************************************
  // initialization of the samples:
  for(t=0;t<numtemp;t++) // traverse the chains
    {
      int first=1, trials=0;

      do // make a new parameter sample without
	// using the data
	{
	  pars[t].numparam=numpar;

	  // cleanup:
	  if(pars[t].param)
	    delete [] pars[t].param;
	  pars[t].param=new double[numpar];

	  if(!startpar || t!=0 || !first)
	    {
	      for(i=0;i<numpar;i++) // traverse the parameters
		pars[t].param[i]=init_par(i);
	    }
	  else
	    {
	      for(i=0;i<numpar;i++)
		pars[t].param[i]=startpar[i];
	    }
	  
	  // update the logprob:
	  log_prob[t]=logprob(pars[t],T[t]);
	  
	  trials++;
	  first=0;
	  if(!silent)
#ifdef MAIN
	    cout << trials << " trials. logprob=" << log_prob[t] << " temp=" <<
	      T[t] << endl;
#else
	  Rcout << trials << " trials. logprob=" << log_prob[t] << " temp=" <<
	    T[t] << std::endl;
#endif // MAIN
	  
	  // check if this logprob makes sense, if not, resample
	} while(!(log_prob[t]>-1e+198 && log_prob[t]<1e+198));
    }
  if(!silent)
#ifdef MAIN
    cout << "done getting initial values" << endl;
#else
  Rcout << "done getting initial values" << std::endl;
#endif // MAIN
  
#ifdef DETAILED_TIMERS
  timers[20][1]=clock();
  timers[20][2]+=(timers[20][1]-timers[20][0]);
#endif // DETAILED_TIMERS

  
  double ***XX=NULL;
  if(dosmooth && !do_realization)
    {
      XX=new double**[num_states];
      for(i=0;i<num_states;i++)
	{
	  XX[i]=new double*[meas_smooth_len];
	  for(k=0;k<meas_smooth_len;k++)
	    {
	      XX[i][k]=new double[num_smooth*numsamples];
	      for(j=0;j<num_smooth*numsamples;j++)
		XX[i][k][j]=0.0;
	    }
	}
    }
  
  
#ifdef DETAILED_TIMERS
  timers[21][0]=clock();
#endif // DETAILED_TIMERS
  
  
  // **************************************************************
  // burn-in phase:
  for(j=0;j<2;j++) // two sub-phases, with and without 
    // adaptation, done twice
    {
      // sample without adaptation:
      for(i=1;i<=burnin/4;i++)
	{
	  if(!silent && talkative_burnin)
#ifdef MAIN
	    cout << j << " 1 " << i << endl;
#else
	  Rcout << j << " 1 " << i << std::endl;
#endif // MAIN
	  newsample(pars, log_prob, &acc, &rw,
		    T, numtemp, swaps);
	}
      
      // update the acceptances:
      for(k=0;k<numpar;k++)
	totalacc.param[k]=0.0;
      
      // sample with adaptation:
      for(i=1;i<=burnin/4;i++)
	{
	  if(!silent && talkative_burnin)
#ifdef MAIN
	    cout << j << " 2 " << i << endl;
#else
	  Rcout << j << " 2 " << i << std::endl;
#endif // MAIN
	    
	  newsample(pars, log_prob, &acc, &rw,
		    T, numtemp, swaps);

	  // update the total amount of acceptances for each parameter:
	  for(k=0;k<numpar;k++)
	    totalacc.param[k]+=acc.param[k];

	  // Adapt according to the acceptance rate each 100th sample:
	  if(i%100==0)
	    {
	      // traverse the aprameters:
	      for(k=0;k<numpar;k++)
		{
		  // if not binary, adapt the standard deviation
		  // so that it goes towards 0.3333:
		  if(par_trans_type[k]!=T_BINARY)
		    rw.param[k] *= exp((totalacc.param[k]/100.0-0.3333)*2.0);

		  // initialize the total amount of accpetances again:
		  totalacc.param[k]=0.0;
		}
	      
	      if(!silent && talkative_burnin)
		{
#ifdef MAIN
		  cout << "RW:";
		  for(k=0;k<numpar;k++)
		    cout << rw.param[k] << " ";
		  cout << endl;
#else
		  Rcout << "RW:";
		  for(k=0;k<numpar;k++)
		    Rcout << rw.param[k] << " ";
		  Rcout << std::endl;
#endif // MAIN
		}
	    }
	}
    }
  
#ifdef DETAILED_TIMERS
  timers[21][1]=clock();
  timers[21][2]+=(timers[21][1]-timers[21][0]);
#endif // DETAILED_TIMERS
  
  
  
  
#ifdef DETAILED_TIMERS
  timers[22][0]=clock();
#endif // DETAILED_TIMERS
  
  
  // **************************************************************
  // MCMC sampling:
  for(i=0;i<numsamples;i++) // traverse the wanted number of samples
    {
      // For each wanted sample, do MCMC iterations 'indep' number
      // of times:
      for(j=0;j<indep;j++)
	newsample(pars, log_prob, &acc, &rw,
		  T, numtemp, swaps,
		  j==(indep-1) ? dosmooth : 0, j==(indep-1) ? do_realization : 0);
      // copy the result to the return array:
      ret[i].copy(pars);
      ll[i]=pars[0].log_lik;
      lp[i]=pars[0].log_prob;

      if(num_series_feed>0)
	{
	  ret[i].iscomplex=is_complex;
	  if(ret[i].iscomplex)
	    {
	      for(j=0;j<10;j++)
		ret[i].cycles[j]=cycle_time[j];
	    }
	  else
	    {
	      for(j=0;j<10;j++)
		ret[i].cycles[j]=MISSING_VALUE;
	    }
	}
      else
	{
	  ret[i].iscomplex=false;
	  for(j=0;j<10;j++)
	    ret[i].cycles[j]=MISSING_VALUE;
	}
      
      
      if(dosmooth && !do_realization)
	{
	  for(k=0;k<meas_smooth_len;k++)
	    {
	      double *xx=new double[num_states], 
		**ss=Make_matrix(num_states,num_states);
	      unsigned int ii,jj;
	      
	      for(ii=0;ii<num_states;ii++)
		xx[ii]=x_k_s_kept[k][ii];
	      for(ii=0;ii<num_states;ii++)
		for(jj=0;jj<num_states;jj++)
		  ss[ii][jj]=P_k_s_kept[k][ii][jj];
		  
	      //check for singularity
	      //double diag=ss[0][0];
	      //double offdiag=ss[0][1];
	      int is_singular=0;
	      /*
		int is_singular=1;
		for(ii=0;ii<num_states && is_singular;ii++)
		for(jj=0;jj<num_states && is_singular;jj++)
		if((ii==jj && !almost_equal(ss[ii][jj], diag)) ||
		(ii!=jj && !almost_equal(ss[ii][jj], offdiag)))
		is_singular=0;
	      */

	      /*
		if(t_k_smooth[k]> -10.045 && t_k_smooth[k]< -10.035)
		cout << is_singular << std::endl;
	      */

	      //if(ml_started)
	      //Rcout << "det ss started" << std::endl;
	      double det=matrix_determinant(ss,num_states);
	      //if(ml_started)
	      //Rcout << "det ss ended" << std::endl;
	      double diagprod=1.0;
	      for(ii=0;ii<num_states;ii++)
		diagprod*=ss[ii][ii];
	      is_singular=(det <= ABSVAL((diagprod*1e-8)) ? 1 : 0);

	      /*
		if(t_k_smooth[k]> -10.045 && t_k_smooth[k]< -10.035)
		{
		cout << is_singular << endl;
		cout << "det=" << det << " " << diagprod << " " << 
		ABSVAL((diagprod*1e-8)) << " " <<
		(det <= ABSVAL((diagprod*1e-8)) ? 1 : 0) << endl;
		cout << "1 " << xx[0] << " " << xx[1] << endl;
		}
	      */

	      // sample:
	      double **new_x;
	      if(!is_singular || num_states==1)
		//new_x=sample_from_multinormal(num_smooth,xx,ss,
		//			      num_states,rptr);
		new_x=multinormal_sample(num_smooth,xx,ss,
					 num_states);
	      else
		{
		  new_x=new double*[num_smooth];
		  for(ii=0;ii<num_smooth;ii++)
		    {
		      new_x[ii]=new double[num_states];
		      double new_x_diag=sqrt(ss[0][0])*get_random_gauss();
		      
		      for(jj=0;jj<num_states;jj++)
			{
			  new_x[ii][jj]=xx[jj]+new_x_diag;
			}
		    }
		}

	      for(j=0;j<num_smooth;j++)
		{
		  for(ii=0;ii<num_states;ii++)
		    XX[ii][k][i*num_smooth+j]=new_x[j][ii];
		}

	      doubledelete(new_x,num_smooth);
	      doubledelete(ss,numsites);
	      delete [] xx;
	    }


	  cleanup_x_and_P(meas_smooth_len);
	}


#ifdef MAIN
      // Show debug info, if wanted:
      if(i%10==0 && !silent)
	{
	  printf("i=%d: ", i);
	  if(num_series_feed>0)
	    printf("%s ", is_complex ? "complex eigenvalues" : "   real eigenvalues");
	  for(j=0;j<numpar;j++)
	    printf("%s=%f ", par_name[j], 
		   invtransform_parameter(pars->param[j], par_trans_type[j]));
	  printf(" l=%g\n", loglik(pars->param));
	  printf("\n");
	}
#endif // MAIN

      if(do_realization)
	{
	  if(!x_k_realized)
	    {
#ifdef MAIN
	      cout << "Skipped realization " << i+1 << endl;
#else
	      Rcout << "Skipped realization " << i+1 << std::endl;
#endif // MAIN
	    }
	  else if(realization_file_start==NULL)
	    {
	      unsigned int i_real=numit_realizations_made;
	      if(i_real<numit_realization)
		{
		  if(!x_k_realized_all)
		    x_k_realized_all=new double**[numit_realization];
		  
		  x_k_realized_all[i_real]=new double*[meas_smooth_len];
		  for(k=0;k<meas_smooth_len;k++)
		    {
		      x_k_realized_all[i_real][k]=new double[num_states];
		      for(j=0;j<num_states;j++)
			x_k_realized_all[i_real][k][j]=x_k_realized[k][j];
		    }
		  numit_realizations_made++;
		}
	    }
	  else
	    {
#ifdef MAIN
	      char filename[1000], cmd[2000];
	      snprintf(filename, 999,
		       "%s_%08d.txt", realization_file_start, i+1);
	      FILE *f=fopen(filename,"w");
	      
	      for(k=0;k<meas_smooth_len;k++)
		{
		  fprintf(f,"%10.4lf ", meas_smooth[k].tm);
		  
		  for(j=0;j<num_states;j++)
		    fprintf(f,"%10.4lf ", x_k_realized[k][j]);
		  fprintf(f,"\n");
		}
	      
	      fclose(f);
	      
	      snprintf(cmd,1500,"gzip %s", filename);
	      int sint=system(cmd);
	      if(sint)
		{
		  printf("system(%s) failed!\n",cmd);
		}
#endif // MAIN
	    }  
	}
    }
  
#ifdef DETAILED_TIMERS
  timers[22][1]=clock();
  timers[22][2]+=(timers[22][1]-timers[22][0]);
#endif // DETAILED_TIMERS
  
  // Show tempering swapping info, if wanted:
  if(!silent)
    {
      // Tempering output:
      if(numtemp>1)
	{
#ifdef MAIN
	  printf("Swaps:\n");
	  for(i=0;i<(numtemp-1);i++)
	    printf("%d<->%d: %d\n", i, i+1, swaps[i]);
#else
	  Rcout << "Swaps:" << std::endl;
	  for(t=0;t<(numtemp-1);t++)
	    Rcout << t << "<->" << t+1 << ": " << swaps[t] << std::endl;
#endif // MAIN
	}
    }


  if(do_importance)
    {
#ifdef DETAILED_TIMERS
  timers[23][0]=clock();
#endif // DETAILED_TIMERS
  
      //Rcout << "loglik importance" << std::endl;

      // **************************************************************
      // Importance sampling for finding the marginal data probability
      // Get mean and correlation for the coefficients 
      // in the different models:
      
      // csamples: will contain the MCMC samples as a 2D array:
      double **csamples=new double*[numsamples];
      // number of non-indicator parameters:
      unsigned int numpar2=(unsigned int)(numpar-numsites*use_indicator);
      double maxprob=-1e+200;

      // Calculate the Bayesian probability estimates for the
      // posterior probability that the indicator=1:
      double prob_indicator[numsites];
      if(use_indicator)
	{
	  for(j=0;j<numsites;j++)
	    {
	      prob_indicator[j]=1.0;
	      for(i=0;i<numsamples;i++)
		prob_indicator[j]+=ret[i].param[j];
	      // Bayesian estimator:
	      prob_indicator[j]/=double(numsamples+2); 
	    }
	}
      
      // Fill the MCMC sample array:
      for(i=0;i<numsamples;i++)
	csamples[i]=ret[i].param+numsites*use_indicator;
      // Fetch the moments:
      double *mu_coefs=mean_of_vectors(numsamples, csamples, numpar2);
      double **sigma_coefs=estimated_variance(numsamples, csamples,
					      mu_coefs, numpar2);
      delete [] csamples;
      
      
      if(!silent)
	{
	  Complex det_eigen=1.0;
	  Complex *lambda=Complex_eigenvalues(sigma_coefs, numpar2);
	  for(i=0;i<numpar;i++)
	    {
	      //cout << "lambda" << i+1 << "=" << lambda[i] << endl;
	      det_eigen*=lambda[i];
	    }
	  double det_orig=matrix_determinant(sigma_coefs, numpar2);
	  Complex logdet=log_matrix_determinant(sigma_coefs, numpar2);
	  Complex detlogc=complex_exp(logdet);
	  double det_log=detlogc.Re();
	  
#ifdef MAIN
	  cout << "det_eigen=" << det_eigen.Re() << "+" << det_eigen.Im() << 
	    "i det_orig=" << det_orig << 
	    " det_log=" << det_log << endl;
#else
	  Rcout << "det_eigen=" << det_eigen.Re() << "+" << det_eigen.Im() << 
	    "i det_orig=" << det_orig << 
	    " det_log=" << det_log << std::endl;
#endif // MAIN

	  delete [] lambda;
	}
      
      // probaiblity sum, contribution and prior sum
      double probsum=0.0;
      double contrib=0.0;
      
      // find standard contrib:
      params par0(numpar);
      if(use_indicator) // if indicator parameters
	{
	  // fetch indicator proposals
	  // and update the proposal density:
	  for(j=0;j<numsites;j++)
	    {
	      if(0.5<prob_indicator[j])
		par0.param[j] = 1.0;
	      else
		par0.param[j] = 0.0;
	    }
	}
      // put the proposals into the proposal parameter structure
      for(j=0;j<numpar2;j++)
	par0.param[j+numsites*use_indicator]=mu_coefs[j];
      double lp0=logprob(par0,1.0,0,0,!silent);
      if(!silent)
#ifdef MAIN
	cout << "lp0=" << lp0 << endl;
#else
      Rcout << "lp0=" << lp0 << std::endl;
#endif // MAIN
      

      unsigned int num_imp=numtemp*(numsamples*indep+burnin)*numpar;
      
      // Importance sampling in order to calculate
      // the model probabilities:
      for(i=0;i<num_imp;i++)
	{
	  // sample from a multinormal proposal distribution having 
	  // the same moments as the posterior MCMC samples:
	  // old
	  /* double **csample=sample_from_multinormal(1, mu_coefs, 
	     sigma_coefs, 
	     numpar2, rptr); */
	  //new 
	  double **csample=multinormal_sample(1, mu_coefs, 
	  				      sigma_coefs, 
	  				      numpar2);

	  // calculate the proposal distribution density:
	  double log_prop_g=pdf_multinormal(*csample, mu_coefs, 
					    sigma_coefs, numpar2, 
					    true);
	  while(!(log_prop_g>=-1e+200 && log_prop_g<=1e+200))
	    {
	      doubledelete(csample,1);
	      csample=multinormal_sample(1, mu_coefs, 
					 sigma_coefs, 
					 numpar2);
	      log_prop_g=pdf_multinormal(*csample, mu_coefs, 
					 sigma_coefs, numpar2, 
					 true);
	    }

	  // make a proposal parameter structure:
	  params par(numpar);
	  if(use_indicator) // if indicator parameters
	    {
	      // fetch indicator proposals
	      // and update the proposal density:
	      for(j=0;j<numsites;j++)
		{
		  if(drand()<prob_indicator[j])
		    {
		      par.param[j] = 1.0;
		      log_prop_g += log(prob_indicator[j]);
		    }
		  else
		    {
		      par.param[j] = 0.0;
		      log_prop_g += log(1.0-prob_indicator[j]);
		    }
		}
	    }
	  
	  // put the proposals into the proposal parameter structure
	  for(j=0;j<numpar2;j++)
	    par.param[j+numsites*use_indicator]=csample[0][j];
	  
	  // fetch the logprob for the proposal:
	  double lp=logprob(par,1.0);
	  
	  double lprob=lp-lp0-log_prop_g;
	  double prob=exp(lprob);
	  
	  if((maxprob>-1e+220 && maxprob<1e+200) && 
	     (maxprob<log_prop_g || !(log_prop_g>=-1e+200 && log_prop_g<1e+200)))
	    maxprob=log_prop_g;
	  
	  // sanity checks, for debugging purposes:
	  if(!(prob<1e+200))
	    {
#ifdef MAIN
	      cout << "noe er galt3!" << endl;
	      cout << "prob=" << prob << " logprob=lp-lp0-log_prop_g=" << lp-lp0-log_prop_g << " lp=" << lp << " lp0=" << lp0 << "log_prop_g=" << log_prop_g << endl;
#else
	      Rcout << "noe er galt3!" << std::endl;
	      Rcout << "prob=" << prob << " logprob=lp-lp0-log_prop_g=" << lp-lp0-log_prop_g << " lp=" << lp << " lp0=" << lp0 << "log_prop_g=" << log_prop_g << std::endl;
#endif // MAIN
	      prob=exp(lp-lp0-log_prop_g);
	    }
	  
	  // update the sum probability and the sum prior:
	  probsum += prob;
	  //prisum += w*prob/expl((long double) ll);
	  
	  // show debug info, if wanted:
	  if(!silent && contrib<prob)
	    {
#ifdef MAIN
	      double ll1=loglik(par.param);
	      double pp=(double) prob;
	      printf("%5d logpropg=%f p=%g lik=%g prior=%g w*p=%g probsum=%g\n", 
		     i, log_prop_g, (double) pp, 
		     exp(ll1), pp/exp(ll1),
		     (double) (prob), 
		     //(double) prisum,
		     (double) probsum);
#endif // MAIN
	      contrib=prob;
	    }
	  
	  // cleanup:
	  doubledelete(csample,1);
	}

      // Monte Carlo estimate for the marginal data density
      probsum/=((double) num_imp);
      // Show the model marginal data density:
#ifdef MAIN
      if(!silent)
	printf("probsum_0=%g\n", double(probsum));
#endif // MAIN
      
      // Show the model marginal data density:
      double lprobsum=log(probsum)+lp0;
      if(!silent)
#ifndef MAIN
        Rcout << "lprobsum=" << lprobsum << std::endl;
#endif // MAIN
      
      if(model_loglik)
	{
	  //(*model_loglik)=log(probsum); 
	  //(*model_loglik)=lp0;
	  (*model_loglik)=lprobsum;
	  //(*model_loglik)=maxprob;
	}
      
      probsum*=exp(lp0);
      if(!silent)
#ifdef MAIN
	printf("probsum=%g\n", double(probsum));
#else
      Rcout << "probsum=" << probsum << std::endl;
#endif // MAIN      
      
      double mean_D=-2.0*find_statistics(ll,numsamples,MEAN);
      double D_mean=-2.0*loglik(mu_coefs);
      double var_D=4.0*find_statistics(ll,numsamples, VARIATION);
      double p_eff1=mean_D-D_mean;
      double p_eff2=0.5*var_D;
      double dic1=mean_D+p_eff1,dic2=mean_D+p_eff2;
      
      if(model_dic1)
	(*model_dic1)=dic1;
      
      if(eff_num_param1)
	(*eff_num_param1)=p_eff1;
      
      if(model_dic2)
	(*model_dic2)=dic2;
      
      if(eff_num_param2)
	(*eff_num_param2)=p_eff2;
      
#ifdef MAIN
      cout << "-DIC1/2 = " << -0.5*dic1 << endl; 
      cout << "( DIC1: " << dic1 << " , mean_D=" << mean_D << " , p_eff1=" << 
	p_eff1 << " , D_mean=" << D_mean << " )" << endl;
      
      cout << "-DIC2/2 = " << -0.5*dic2 << endl; 
      cout << "( DIC2: " << dic2 << " , mean_D=" << mean_D << " , p_eff2=" << 
	p_eff2 << " , var_D=" << var_D << " )" << endl;
#endif // MAIN
      
      delete [] mu_coefs;
      doubledelete(sigma_coefs,numpar2);
      
#ifdef DETAILED_TIMERS
      timers[23][1]=clock();
      timers[23][2]+=(timers[23][1]-timers[23][0]);
#endif // DETAILED_TIMERS
  
    }
  else
    {
      if(model_loglik)
	*model_loglik=MISSING_VALUE;
    }
  
  if(x_)
    *x_=XX;
  
  // return the samples:
  return ret;
}
	

// *************************************************
// show_parameters: Uses the external programs
// 'vvgraph' and 'histogramme' from the hydrasub
// package to show the sampling for a given
// parameter. If these are not present, this
// methods should be set on silent (only spacing
// between independent samples are shown).
// par: parameter samples
// N: number of samples
// parname: name of parameter
// silent: if set, turns on silent modus
// *************************************************
void show_parameter(double *par, int N, char *parname, int silent, char *filestart) 
{
#ifdef MAIN
  char cmd[1000], filename[1000]; // file name string
  int i;
#endif // MAIN
  
  int len=N; 
  // calculate one-step autocorrelation:
  double rho=get_auto_correlation(par, len);
  // number of independent samples, using 
  // AR(1) model and Gelman's formula for independent samples:
  // double n_indep=len/2.0/(0.5+rho/(1.0-rho));
  // spacing=N/n_indep
  double spacing=2.0*(0.5+rho/(1.0-rho));

  // Show spacing:
#ifdef MAIN
  cout << parname << " - spacing between independent samples:" << 
    spacing << endl;
#else
  Rcout << parname << " - spacing between independent samples:" << 
    spacing << std::endl;
#endif // MAIN

#ifdef MAIN
  
  cout << parname << " mean=" << find_statistics(par,N,MEAN);
  cout << " median=" << find_statistics(par,N,MEDIAN);
  cout << " 95% post.cred.int=(" << 
    find_statistics(par,N,PERCENTILE_2_5) << " , " <<
    find_statistics(par,N,PERCENTILE_97_5) << ")" << endl;
  
  // if silent, this is all that is to be done:
  if(silent)
    return;

  FILE *p; // file pointer
  // show sampling history using 'vvgraph':
  if(!filestart || !*filestart)
    {
      snprintf(cmd, 999, "vvgraph -x \"%s\" 2> /dev/null > /dev/null", parname);
      p=popen(cmd,"w");
    }
  else
    {
      snprintf(filename, 999, "%s_%s_ts.txt", filestart, parname);
      p=fopen(filename,"w");
    }
  fprintf(p,"# Column 1: %s\n",  parname);
  fprintf(p,"###################\n");
  for(i=0;i<len;i++)
    fprintf(p,"%d %lf\n", i+1, par[i]);
  if(!filestart || !*filestart)
    pclose(p);
  else
    fclose(p);

  // show histogram using 'histogramme':
  if(!filestart || !*filestart)
    {
      snprintf(cmd, 999,
	       "histogramme -x \"%s\" -t \"%s\" 2> /dev/null > /dev/null", 
	       parname,  parname); 
      p=popen(cmd, "w");
    }
  else
    {
      snprintf(filename, 999, "%s_%s_hist.txt", filestart, parname);
      p=fopen(filename,"w");
    }
  for(i=0;i<len;i++)
    fprintf(p,"%lf\n", par[i]);
  if(!filestart || !*filestart)
    pclose(p);
  else
    fclose(p);
#endif // MAIN
}


#ifdef MAIN
// ******************************************************
// show_scatter: Shows scatterplots of the samples for
// two parameter.
// par1: parameter samples for parameter 1
// par2: parameter samples for parameter 2
// N: number of samples
// parname1: parameter name 1
// parname2: parameter name 2
// ******************************************************
void show_scatter(double *par1, double *par2, int N, 
		  char *parname1, char *parname2, char *filestart)
{
  FILE *p; // file pointer
  char cmd[1000], filename[1000]; // command string + file name string
  int i, len=N;
  
  // show scatterplot using the 'vvgraph' program:

  if(!filestart || !*filestart)
    {
      snprintf(cmd, 999, "vvgraph -x \"%s\" -y \"%s\"", parname1, parname2);
      p=popen(cmd,"w");
    }
  else
    {
      snprintf(filename, 999,
	       "%s_%s_%s_scatter.txt", filestart, parname1,parname2);
      p=fopen(filename,"w");
    }
  fprintf(p, "# Column 1: %s vs %s\n", parname1, parname2);
  fprintf(p, "# Column 1 - type: dot\n");
  fprintf(p, "#####################\n");
  for(i=0;i<len;i++)
    fprintf(p,"%lf %lf\n", par1[i], par2[i]);
  if(!filestart || !*filestart)
    pclose(p);
  else
    fclose(p);
}
#endif // MAIN


#ifndef MAIN // Only for when compiling as R package

#include <Rcpp.h>
#include <RcppCommon.h>



const double loglikwrapper(const NumericVector &vals)
{
  int len=vals.size();
  double *val2=new double[len];

  for(int i=0;i<len;i++)
    val2[i]=vals[i];

  double ret=-minusloglik(val2);
  delete [] val2;

  return -ret;
}


double *partial_par=NULL;
const double partial_loglikwrapper(const NumericVector &vals)
{
  if(partial_par==NULL)
    return MISSING_VALUE;
  
  double *val2=new double[numpar];

  unsigned int i=0,j;
  for(j=0;j<numpar;j++)
    {
      if(partial_par[j]!=MISSING_VALUE)
	val2[j]=partial_par[j];
      else
	{
	  val2[j]=vals[i];
	  i++;
	}
    }
  
  double ret=-minusloglik(val2);
  delete [] val2;
  
  return -ret;
}


RcppExport SEXP layeranalyzer(SEXP input,SEXP num_MCMC ,SEXP Burnin,
			      SEXP Spacing,SEXP NumTemp,
			      SEXP do_model_likelihood,
			      SEXP do_maximum_likelihood,
			      SEXP maximum_likelihood_numstart,
			      SEXP SilentMode, SEXP TalkativeBurnin,
			      SEXP TalkativeLikelihood,
			      SEXP IdStrategy, SEXP UseStationarySdev,
			      SEXP TempGround, SEXP UseHalfLives,
			      SEXP ReturnMCMC, 
                              SEXP causal, SEXP causal_symmetric, SEXP corr,
			      SEXP smooth_specs, SEXP realization_specs,
			      SEXP ReturnResiduals,
			      SEXP input_param_values)
{  
  reset_global_variables();
  ser=new series[100];  
  id_strategy=ID_NONE; // ID_SUB_LOWER;
  
  GetRNGstate();
  int static first=1;
  if(first)
    randify();
  first=0;
  
  int num_mcmc=as<int>(num_MCMC);
  int burnin=as<int>(Burnin);
  int spacing=as<int>(Spacing);
  unsigned int numtemp=as<unsigned int>(NumTemp);
  int do_importance=as<int>(do_model_likelihood);
  int do_ml=as<int>(do_maximum_likelihood);
  if(do_ml)
    do_importance=0;
  int num_optim=as<int>(maximum_likelihood_numstart);
  silent=as<int>(SilentMode);
  talkative_burnin=as<int>(TalkativeBurnin);
  talkative_likelihood=as<int>(TalkativeLikelihood);
  id_strategy=(INDENTIFICATION_PRIOR_HANDLING)as<int>(IdStrategy);
  int use_stationary_sdev=as<int>(UseStationarySdev);
  double T_ground=as<double>(TempGround);
  int use_half_times=as<int>(UseHalfLives);
  int return_mcmc=as<int>(ReturnMCMC);
  int return_residuals=as<int>(ReturnResiduals);
  
  NumericMatrix in_pars=as<NumericMatrix>(input_param_values); 
  int inpars_numsets=in_pars.nrow();
  int inpars_numpars=in_pars.ncol();
  bool no_inpars = (inpars_numsets==1 && inpars_numpars==1) ? true : false;
  
  double **resids=NULL, *resids_time, **prior_expected_values=NULL;
  int resid_numcol=0, resid_len=0;
  
  // DEBUG: Rcout << "return_residuals=" << return_residuals << std::endl;
  
  List smoothspecs=as<List>(smooth_specs);
  int dosmooth=as<int>(smoothspecs["do.smoothing"]);
  double smooth_diff=MISSING_VALUE; 
  int num_smooth_per_mcmc=10;
  int do_return_smoothing_samples=0;
  int do_start_end=0, do_start_end_dt=0; 
  double smooth_start=MISSING_VALUE,smooth_end=MISSING_VALUE;
  HydDateTime smooth_start_dt=NoHydDateTime, smooth_end_dt=NoHydDateTime;
  if(dosmooth)
    {
      smooth_diff=as<double>(smoothspecs["smoothing.time.diff"]);
      do_start_end=as<int>(smoothspecs["start.end.given"]);
      do_return_smoothing_samples=as<int>(smoothspecs["do.return.smoothing.samples"]);
      num_smooth_per_mcmc=as<int>(smoothspecs["num.smooth.per.mcmc"]);
      if(do_start_end)
	{
	  do_start_end_dt=as<int>(smoothspecs["start.end.datetime"]);
	  if(do_start_end_dt)
	    {
	      HydDateTime dt1(1970,1,1,0,0);
	      dt1+=(int) as<int>(smoothspecs["smoothing.start.dt"]);
	      smooth_start_dt=dt1;
	      
	      HydDateTime dt2(1970,1,1,0,0);
	      dt2+=(int) as<int>(smoothspecs["smoothing.end.dt"]);
	      smooth_end_dt=dt2;
	    }
	  else
	    {
	      smooth_start=as<double>(smoothspecs["smoothing.start"]);
	      smooth_end=as<double>(smoothspecs["smoothing.end"]); 
	    }
	}
      
      num_smooth=num_smooth_per_mcmc;
      dosmooth=1;
    }
  
  List realizespecs=as<List>(realization_specs);
  int do_realizations=as<int>(realizespecs["do.realizations"]);
  int do_realize_start_end=0; // , do_realize_start_end_dt=0; 
  double realize_start=MISSING_VALUE,realize_end=MISSING_VALUE;
  HydDateTime realize_start_dt=NoHydDateTime, realize_end_dt=NoHydDateTime;
  double realize_diff=MISSING_VALUE; 
  real_strat=NO_CENSORING;
  numit_realization=0;
  if(do_realizations)
    {
      numit_realization=as<int>(realizespecs["num.realizations"]);
      
      std::string rstrat_string=as<std::string>(realizespecs["strategy"]);
      char *rstrat_char=(char *) rstrat_string.c_str();
      
      if(rstrat_char[0]=='A' || rstrat_char[0]=='a') // ascending restriction
	real_strat=ASCENDING_CENSORING;
      else if(rstrat_char[0]=='D' || rstrat_char[0]=='d') // descending restriction
	real_strat=DESCENDING_CENSORING;
      
      realize_diff=as<double>(realizespecs["realization.time.diff"]);
      do_realize_start_end=as<int>(realizespecs["start.end.given"]);
      if(do_realize_start_end)
	{
	  //do_realize_start_end_dt=as<int>(realizespecs["start.end.datetime"]);
	  if(do_start_end_dt)
	    {
	      HydDateTime dt1(1970,1,1,0,0);
	      dt1+=(int) as<int>(realizespecs["realization.start.dt"]);
	      realize_start_dt=dt1;
	      
	      HydDateTime dt2(1970,1,1,0,0);
	      dt2+=(int) as<int>(realizespecs["realization.end.dt"]);
	      realize_end_dt=dt2;
	    }
	  else
	    {
	      realize_start=as<double>(realizespecs["realization.start"]);
	      realize_end=as<double>(realizespecs["realization.end"]); 
	    }
	}
      
    }

  
  List inall(input);
  num_series=inall.size();
  
  if(!silent)
    Rcout << inall.size() << std::endl;
  
  unsigned int s,i,l;
  int j;
  for(s=0;s<num_series;s++)
    {
      List inlist=as<List>(inall[s]);
      
      List indata=as<List>(inlist["timeseries"]);
      
      NumericVector loglik(1);
      
      loglik=MISSING_VALUE;
      
      // int N=inlist.size();
      //NumericMatrix X=as<NumericMatrix>(inlist["datamatrix"]); 
      
      NumericVector X_time=as<NumericVector>(indata["time"]);
      unsigned int len_time=X_time.size();
      
      NumericVector X_value=as<NumericVector>(indata["value"]);
      unsigned int len_value=X_value.size();
      
      if(len_time!=len_value)
	{
	  Rcout << "Length of 'time' array ("<< len_time <<
	    ")!= length of 'value' array ("
		<< len_value << ")!" << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      unsigned int len=len_time;
      
      std::string name=as<std::string>(indata["name"]);
      char *seriename=(char *) name.c_str();
      
      int num_layers=as<int>(inlist["numlayers"]);
      int nopull=as<int>(inlist["no.pull"]); 
      int lintime=as<int>(inlist["lin.time"]); 
      IntegerVector timeint_layers=as<IntegerVector>(inlist["time.integral"]);
      unsigned int num_timeint=timeint_layers.size();
      int regmu=as<int>(inlist["regional.mu"]); 
      int reglint=as<int>(inlist["regional.lin.time"]); 
      IntegerVector regpull_layers=as<IntegerVector>(inlist["regional.pull"]);
      unsigned int num_regpull=regpull_layers.size();
      IntegerVector regsigma_layers=as<IntegerVector>(inlist["regional.sigma"]);
      unsigned int num_regsigma=regsigma_layers.size();
      IntegerVector corrsigma_layers=as<IntegerVector>(inlist["correlated.sigma"]);
      unsigned int num_corrsigma=corrsigma_layers.size();
      IntegerVector pcorrsigma_layers=as<IntegerVector>(inlist["pairwise.correlated.sigma"]);
      unsigned int num_pcorrsigma=pcorrsigma_layers.size();
      IntegerVector nosigma_layers=as<IntegerVector>(inlist["no.sigma"]);
      unsigned int num_nosigma=nosigma_layers.size();
      IntegerVector onedimsigma_layers=as<IntegerVector>(inlist["one.dim.sigma"]);
      unsigned int num_onedimsigma=onedimsigma_layers.size();
      IntegerVector groupingsigma_layers=as<IntegerVector>(inlist["grouping.sigma"]);
      unsigned int num_groupingsigma=groupingsigma_layers.size();
      IntegerVector removesigma_layers=as<IntegerVector>(inlist["remove.sigma"]);
      unsigned int num_removesigma=removesigma_layers.size();
      IntegerVector diffsigma_layers=as<IntegerVector>(inlist["differentiate.sigma"]);
      unsigned int num_diffsigma=diffsigma_layers.size();
      IntegerVector diffpull_layers=as<IntegerVector>(inlist["differentiate.pull"]);
      unsigned int num_diffpull=diffpull_layers.size();
      int diffmu=as<int>(inlist["differentiate.mu"]); 
      int difflint=as<int>(inlist["differentiate.lin.time"]); 
      int init0=as<int>(inlist["init.0"]); 
      double init_time=as<double>(inlist["init.time"]); 
      int init_datetime=as<int>(inlist["init.datetime"]);
      int init_same_sites=as<int>(inlist["init.same.sites"]); 
      int init_same_layers=as<int>(inlist["init.same.layers"]); 
      NumericVector init_specified=as<NumericVector>(inlist["init.specified"]);
      unsigned int num_init_specified=init_specified.size();
      int allow_pos_pull=as<int>(inlist["allow.pos.pull"]);
      NumericVector periods=as<NumericVector>(inlist["period"]);
      unsigned int num_periods=periods.size();
      int is_datetime=as<int>(indata["is.datetime"]);
      
      unsigned int n=len;
      
      NumericVector sd=as<NumericVector>(indata["std.dev"]);
      unsigned int len_sd=sd.size();
      double *sd_meas=NULL;
      if(len_sd==len)
	{
	  sd_meas=new double[len];
	  for(i=0;i<len;i++)
	    sd_meas[i]=sd[i];
	}
      
      IntegerVector num_meas_per_value=as<IntegerVector>(indata["num.meas.per.value"]);
      unsigned int len_n=num_meas_per_value.size();
      int *n_meas=NULL;
      if(len_n==len)
	{
	  n_meas=new int[len];
	  for(i=0;i<len;i++)
	    n_meas[i]=num_meas_per_value[i];
	}
      
      numsites=1;
      IntegerVector site=as<IntegerVector>(indata["site"]);
      unsigned int len_sites=site.size();
      int *site_meas=NULL;
      if(len_sites==len)
	{
	  site_meas=new int[len];
	  for(i=0;i<len;i++)
	    {
	      site_meas[i]=site[i];
	      if((site_meas[i]+1)>(int)numsites)
		numsites=(unsigned int) (site_meas[i]+1);
	    }
	}  
      
      double *Xtime=new double[n], *Xval=new double[n], meanval=0.0;
      j=0;
      for(i=0;i<len;i++)
	{
	  Xtime[i]=X_time[i];
	  Xval[i]=X_value[i];
	  if(X_value[i]!=MISSING_VALUE)
	    {
	      meanval+=X_value[i];
	      j++;
	    }
	}
      meanval/=double(j);
      
      
      List currprior=as<List>(inlist["prior"]);
      NumericVector prior_mu=as<NumericVector>(currprior["mu"]);
      NumericVector prior_dt=as<NumericVector>(currprior["dt"]);
      NumericVector prior_s=as<NumericVector>(currprior["s"]);
      NumericVector prior_lin=as<NumericVector>(currprior["lin"]);
      NumericVector prior_beta=as<NumericVector>(currprior["beta"]);
      NumericVector prior_init=as<NumericVector>(currprior["init"]);
      NumericVector prior_obs=as<NumericVector>(currprior["obs"]);
      int prior_islog=as<int>(currprior["islog"]);
      if(!silent)
	Rcout << "prior_obs1=" << prior_obs[0] << " prior_obs2=" << prior_obs[1] << std::endl;
      prior *newprior=new prior(prior_islog, prior_mu[0],prior_mu[1],
				prior_dt[0],prior_dt[1],
				prior_s[0],prior_s[1],
				prior_lin[0],prior_lin[1],
				prior_beta[0],prior_beta[1],
				prior_init[0],prior_init[1],
				prior_obs[0],prior_obs[1],
				meanval);
      if(!silent)
	{
	  Rcout << "prior_obs1=" << newprior->os_1 << 
	    " prior_obs2=" << newprior->os_2 << std::endl;
	  
	  Rcout << "num_layers=" << num_layers << std::endl;
	}
      
      
      ser[s].set_series(seriename,s,
			Xtime,Xval,n,num_layers,newprior,is_datetime ? true : false,
			sd_meas, n_meas, site_meas);
      
      delete [] Xtime;
      delete [] Xval;
      delete newprior;
      if(sd_meas)
	delete [] sd_meas;
      if(n_meas)
	delete [] n_meas;
      if(site_meas)
	delete [] site_meas;
      
      if(nopull)
	ser[s].no_pull_lower=1;
      if(lintime)
	ser[s].linear_time_dep=true;
      for(i=0;i<num_timeint;i++)
	{
	  //cout << "timeint" << i+1 << ":" << timeint_layers[i] << endl;
	  if(timeint_layers[i])
	    {
	      if(timeint_layers[i]>(num_layers-1))
		{
		  Rcout << "Can't put time integration on the "
		    "lowest layer (or lower)" << std::endl;
		  return NULL;
		}
	      ser[s].time_integral[timeint_layers[i]-1]=1;
	      id_strategy=ID_NONE;
	      ser[s].no_sigma[timeint_layers[i]-1]=1;
	    }
	}
      if(regmu)
	ser[s].regional_mu=1;
      if(reglint)
	{
	  ser[s].regional_lin_t=1;
	  ser[s].linear_time_dep=true;
	}
      for(i=0;i<num_regpull;i++)
	if(regpull_layers[i])
	  ser[s].regional_pull[regpull_layers[i]-1]=1;
      for(i=0;i<num_regsigma;i++)
	if(regsigma_layers[i])
	  ser[s].regional_sigma[regsigma_layers[i]-1]=1;
      for(i=0;i<num_corrsigma;i++)
	if(corrsigma_layers[i])
	  ser[s].sigma_correlated[corrsigma_layers[i]-1]=1;
      for(i=0;i<num_pcorrsigma;i++)
	if(pcorrsigma_layers[i])
	  {
	    ser[s].sigma_pairwise_correlated[pcorrsigma_layers[i]-1]=1;
	    some_pairwise_correlations=true;
	  }
      for(i=0;i<num_nosigma;i++)
	if(nosigma_layers[i])
	  ser[s].no_sigma[nosigma_layers[i]-1]=1;
      for(i=0;i<num_onedimsigma;i++)
	if(onedimsigma_layers[i])
	  ser[s].sigma_1dim[onedimsigma_layers[i]-1]=1;
      for(i=0;i<num_groupingsigma;i++)
	if(groupingsigma_layers[i])
	  ser[s].indicator_corr[groupingsigma_layers[i]-1]=1;
      for(i=0;i<num_removesigma;i++)
	if(removesigma_layers[i])
	  ser[s].indicator_corr2[removesigma_layers[i]-1]=1;
      for(i=0;i<num_diffsigma;i++)
	if(diffsigma_layers[i])
	  ser[s].indicator_sigma[diffsigma_layers[i]-1]=1;
      for(i=0;i<num_diffpull;i++)
	if(diffpull_layers[i])
	  ser[s].indicator_pull[diffpull_layers[i]-1]=1;
      if(diffmu)
	ser[s].indicator_mu=1;
      if(difflint)
	{
	  ser[s].indicator_lin_t=1;
	  ser[s].linear_time_dep=true;
	}
      if(init0)
	ser[s].init_treatment=1;
      if(init_time!=MISSING_VALUE)
	{
	  ser[s].init_treatment=1;
	  ser[s].init_time=init_time;
	  if(init_datetime)
	    {
	      HydDateTime dt(1970,1,1,0,0);
	      dt+=(int) init_time;
	      ser[s].init_dt=dt;
	    }
	}
      if(init_same_sites)
	ser[s].regional_init=0;
      else
	ser[s].regional_init=1;
      if(init_same_layers)
	ser[s].layered_init=0;
      else
	ser[s].layered_init=1;
      if(num_init_specified!=2)
	Rcout << "Input vector for init.specified should be of size 2!" << std::endl;
      else if(init_specified[0]!=MISSING_VALUE)
	{
	  ser[s].init_treatment=1;
	  ser[s].init_time=init_specified[0];
	  HydDateTime dt(1970,1,1,0,0);
	  dt+=(int) ser[s].init_time;
	  ser[s].init_dt=dt; 
	  if(init_specified[1]!=MISSING_VALUE)
	    ser[s].init_value=init_specified[1];
	}
      if(allow_pos_pull)
	{
	  ser[s].allow_positive_pulls=true;
	  id_strategy=ID_NONE;
	}
      for(i=0;i<num_periods;i++)
	if(periods[i]!=MISSING_VALUE && periods[i]>0.0)
	  {
	    ser[s].Tper[ser[s].num_per]=periods[i];
	    ser[s].num_per++;
	  }
    }
  
  num_states=0;
  num_tot_layers=0;
  for(s=0;s<num_series;s++)
    {
      for(l=0;l<ser[s].numlayers;l++)
	for(j=0;j<(int)numsites;j++)
	  {
	    state_series[num_states+l*numsites+j]=s;
	    state_layer[num_states+l*numsites+j]=l;
	  }
      
      series_state_start[s]=num_states;
      num_states += ser[s].numlayers*numsites;
      num_tot_layers+=ser[s].numlayers;
    }  
  
  meas_tot_len=0;
  for(s=0;s<num_series;s++)
    {
      meas_tot_len += ser[s].meas_len;
      for(i=0;i<ser[s].meas_len;i++)
	ser[s].meas[i].index=series_state_start[s]+ser[s].meas[i].site;
    }
  
    
  series_measurements *meas_buffer=new series_measurements[meas_tot_len];	  
  j=0;
  for(s=0;s<num_series;s++)
    for(i=0;i<ser[s].meas_len;i++)
      meas_buffer[j++]=ser[s].meas[i];
  qsort(meas_buffer,(size_t) meas_tot_len,sizeof(series_measurements), 
	compare_meas);
  
  meas_tot=new measurement_cluster[meas_tot_len];
  j=-1;
  for(i=0;i<meas_tot_len;i++)
    {
      if(j<0 || meas_tot[j].tm!=meas_buffer[i].tm)
	j++;
      if(meas_buffer[i].meanval!=MISSING_VALUE)
	meas_tot[j].add_measurement(meas_buffer[i]);
      else if(meas_tot[j].num_measurements==0)
	{
	  meas_tot[j].tm=meas_buffer[i].tm;
	  meas_tot[j].dt=meas_buffer[i].dt;
	}
    }
  meas_tot_len=j+1;
  
  
  
  // Make measurment clusters for smoothed/realized time series:
  if(dosmooth || do_realizations)
    {
      double t_start=dosmooth ? smooth_start : realize_start, 
	t_end=dosmooth ? smooth_end : realize_end;
      //HydDateTime dt_start= smooth_start_dt, dt_end=smooth_end_dt;
      HydDateTime ref(1970);
      double diff=dosmooth ? smooth_diff : realize_diff;
      
      if(diff<=0.0)
	{
	  meas_smooth_len=meas_tot_len;
	  meas_smooth=new measurement_cluster[meas_tot_len];
	  for(i=0;i<meas_tot_len;i++)
	    meas_smooth[i].copy(&(meas_tot[i]));
	}
      else
	{
	  meas_smooth_len=0;
	  double t1=(t_start==MISSING_VALUE ? meas_tot[0].tm-1.0 : t_start);
	  double t2=(t_end==MISSING_VALUE ? 
		     MAXIM(0.0,meas_tot[meas_tot_len-1].tm+1.0) :
		     t_end);
	  double t=t1;
	  
	  if(meas_tot[0].dt!=NoHydDateTime)
	    {
	      if(t_start==MISSING_VALUE)
		t1=meas_tot[0].tm-
		  0.1*double(meas_tot[meas_tot_len-1].dt.
			     difference_minutes(meas_tot[0].dt))*60.0;

	      if(t_end==MISSING_VALUE)
		t2=meas_tot[meas_tot_len-1].tm+
		  0.1*double(meas_tot[meas_tot_len-1].dt.
			     difference_minutes(meas_tot[0].dt))*60.0;
	      t=t1;
	    }
	  
	  for(t=t1;t<=t2;t+=diff)
	    meas_smooth_len++;
	  
	  meas_smooth=new measurement_cluster[meas_smooth_len+meas_tot_len+10];
	  unsigned int j=0,k=0;
	  
	  i=0;
	  t=MINIM(t1,meas_tot[0].tm);
	  while(i<meas_smooth_len+meas_tot_len+1 && 
		(t<MAXIM(t2,meas_tot[meas_tot_len-1].tm)))
	    {
	      while(j<meas_tot_len && meas_tot[j].tm<t)
		{
		  //cout << k << " av " << meas_smooth_len+meas_tot_len+10 << " "
		  //   << j << " av " << meas_tot_len << endl;
		  
		  meas_smooth[k].copy(&(meas_tot[j]));
		  //cout << k << " " << meas_smooth[k].dt << " " << 
		  //meas_tot[j].dt << " " << meas_tot[0].dt << endl;
		  k++;
		  j++;
		}
	      
	      if((j>=meas_tot_len || meas_tot[j].tm>t) &&
		 (t>=t1 && t<=t2))
		{
		  measurement_cluster buff(t);
		  meas_smooth[k].copy(&buff);
		  if(meas_tot[0].dt!=NoHydDateTime)
		    {
		      meas_smooth[k].dt=ref+((long int)t);
		      //cout << ref << " " << t << " " << ((long int) t) << endl;
		      //cout << k << " " << meas_smooth[k].dt << " " << 
		      //meas_tot[0].dt << endl;
		      if(!meas_smooth[k].dt.legal())
			{
			  Rcout << "Illegal date-time:" << ref << " " << 
			    t << " " << ((long int)t) << 
			    " " << meas_tot[0].dt << std::endl;
			  return NULL;
			}
		    }
		  k++;
		}
	      
	      i++;
	      t=t1+diff*double(i);
	    }
	  
	  meas_smooth_len=k;
	}
    }
  
    

  /* DEBUG
     if(return_residuals)
     {
     Rcout << "right after reading:" << endl;
     for(int k=0;k<meas_tot_len;k++)
     meas_tot[k].print();
     Rcout << "later:" << endl;
     }
  */
  
  // Read causal connections:
  feed_from_series=new int[10000];
  feed_to_series=new int[10000];
  feed_from_layer=new int[10000];
  feed_to_layer=new int[10000];
  beta_feed=new double[10000];
  feed_symmetric=new int[10000];
  
  IntegerVector causal_pairs=as<IntegerVector>(causal);
  num_series_feed=causal_pairs.size()/4;
  
  if(!silent)
    {
      Rcout << causal_pairs.size() << std::endl;
      for(i=0;i<causal_pairs.size();i++)
        Rcout << "i=" << i << " f=" << causal_pairs[i] << std::endl;
    }
  
  for(i=0;i<num_series_feed;i++)
    {
      feed_from_series[i]=causal_pairs[4*i]-1;
      if(feed_from_series[i]>=(int)num_series)
	{
	  Rcout << "Series index too high for 'from' series:" << 
	    feed_from_series[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      if(feed_from_series[i]<0)
	{
	  Rcout << "Series index too low for 'from' series:" << 
	    feed_from_series[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      feed_from_layer[i]=causal_pairs[4*i+1]-1;
      if(feed_from_layer[i]>=
	 (int)ser[feed_from_series[i]].numlayers)
	{
	  Rcout << "Layer index too high for 'from' series:" << 
	    feed_from_series[i]+1 << ":" << feed_from_layer[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      if(feed_from_layer[i]<0)
	{
	  Rcout << "Layer index too low for 'from' series:" << 
	    feed_from_series[i]+1 << ":" << feed_from_layer[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      feed_to_series[i]=causal_pairs[4*i+2]-1;
      if(feed_to_series[i]>=(int)num_series)
	{
	  Rcout << "Series index too high for 'to' series:" << 
	    feed_to_series[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      if(feed_to_series[i]<0)
	{
	  Rcout << "Series index too low for 'to' series:" << 
	    feed_to_series[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      feed_to_layer[i]=causal_pairs[4*i+3]-1;
      if(feed_to_layer[i]>=
	 (int)ser[feed_to_series[i]].numlayers)
	{
	  Rcout << "Layer index too high for 'to' series:" << 
	    feed_to_series[i]+1 << ":" << feed_to_layer[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      if(feed_to_layer[i]<0)
	{
	  Rcout << "Layer index too low for 'to' series:" << 
	    feed_to_series[i]+1 << ":" << feed_to_layer[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      feed_symmetric[i]=0;
    }
  
  IntegerVector causal_sym_pairs=as<IntegerVector>(causal_symmetric);
  int numsym=causal_sym_pairs.size()/4;
  for(i=num_series_feed;i<num_series_feed+numsym;i++)
    {
      feed_from_series[i]=causal_sym_pairs[4*(i-num_series_feed)]-1;
      if(feed_from_series[i]>=(int)num_series)
	{
	  Rcout << "Series index too high for 'from' series:" << 
	    feed_from_series[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      if(feed_from_series[i]<0)
	{
	  Rcout << "Series index too low for 'from' series:" << 
	    feed_from_series[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      feed_from_layer[i]=causal_sym_pairs[4*(i-num_series_feed)+1]-1;
      if(feed_from_layer[i]>=
	 (int)ser[feed_from_series[i]].numlayers)
	{
	  Rcout << "Layer index too high for 'from' series:" << 
	    feed_from_series[i]+1 << ":" << feed_from_layer[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      if(feed_from_layer[i]<0)
	{
	  Rcout << "Layer index too low for 'from' series:" << 
	    feed_from_series[i]+1 << ":" << feed_from_layer[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      feed_to_series[i]=causal_sym_pairs[4*(i-num_series_feed)+2]-1;
      if(feed_to_series[i]>=(int)num_series)
	{
	  Rcout << "Series index too high for 'to' series:" << 
	    feed_to_series[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      if(feed_to_series[i]<0)
	{
	  Rcout << "Series index too low for 'to' series:" << 
	    feed_to_series[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      feed_to_layer[i]=causal_sym_pairs[4*(i-num_series_feed)+3]-1;
      if(feed_to_layer[i]>=
	 (int)ser[feed_to_series[i]].numlayers)
	{
	  Rcout << "Layer index too high for 'to' series:" << 
	    feed_to_series[i]+1 << ":" << feed_to_layer[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      if(feed_to_layer[i]<0)
	{
	  Rcout << "Layer index too low for 'to' series:" << 
	    feed_to_series[i]+1 << ":" << feed_to_layer[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      feed_symmetric[i]=1;
    }
  num_series_feed+=numsym;
  
  
  IntegerVector corr_pairs=as<IntegerVector>(corr);
  num_series_corr=corr_pairs.size()/4;
  if(!silent)
    {
      Rcout << corr_pairs.size() << std::endl;
      for(i=0;i<corr_pairs.size();i++)
	Rcout << "i=" << i << " c=" << corr_pairs[i] << std::endl;
    }
  
  corr_from_series=new int[10000];
  corr_to_series=new int[10000];
  corr_from_layer=new int[10000];
  corr_to_layer=new int[10000];
  corr_from_index=new int[10000];
  corr_to_index=new int[10000];
  series_corr=new double[10000];
  
  for(i=0;i<num_series_corr;i++)
    {
      corr_from_series[i]=corr_pairs[4*i]-1;
      if(corr_from_series[i]>=(int)num_series)
	{
	  Rcout << "Series index too high for 'from' series:" << 
	    corr_from_series[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      if(corr_from_series[i]<0)
	{
	  Rcout << "Series index too low for 'from' series:" << 
	    corr_from_series[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      corr_from_layer[i]=corr_pairs[4*i+1]-1;
      if(corr_from_layer[i]>=
	 (int)ser[corr_from_series[i]].numlayers)
	{
	  Rcout << "Layer index too high for 'from' series:" << 
	    corr_from_series[i]+1 << ":" << corr_from_layer[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      if(corr_from_layer[i]<0)
	{
	  Rcout << "Layer index too low for 'from' series:" << 
	    corr_from_series[i]+1 << ":" << corr_from_layer[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      corr_to_series[i]=corr_pairs[4*i+2]-1;
      if(corr_to_series[i]>=(int)num_series)
	{
	  Rcout << "Series index too high for 'to' series:" << 
	    corr_to_series[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      if(corr_to_series[i]<0)
	{
	  Rcout << "Series index too low for 'to' series:" << 
	    corr_to_series[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      corr_to_layer[i]=corr_pairs[4*i+3]-1;
      if(corr_to_layer[i]>=
	 (int)ser[corr_to_series[i]].numlayers)
	{
	  Rcout << "Layer index too high for 'to' series:" << 
	    corr_to_series[i]+1 << ":" << corr_to_layer[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
      if(corr_to_layer[i]<0)
	{
	  Rcout << "Layer index too low for 'to' series:" << 
	    corr_to_series[i]+1 << ":" << corr_to_layer[i]+1 << std::endl;
	  PutRNGstate();
	  return NULL;
	}
    }
  
  if(num_series_corr>0)
    {
      for(i=0;i<num_series_corr;i++)
	{
	  int i1=series_state_start[corr_from_series[i]]+numsites*corr_from_layer[i];
	  int i2=series_state_start[corr_to_series[i]]+numsites*corr_to_layer[i]; 
	  
	  corr_from_index[i]=i1/numsites;
	  corr_to_index[i]=i2/numsites;
	}
    }
  
  if(num_series_corr>1)
    {
      // Estimate the probability of a correlation drawn from the prior
      // resulting in a positively definite covariance matrix.
      // Used for correcting this prior to this restriction.
      p_pos_series_sigma2=p_positive_seriescorr(10000);
    }
  else
    p_pos_series_sigma2=1.0;
  
  
    
  
  
  double ***x=NULL; // smoothing results from the latent processes
  double ***P_k_smooth=NULL, **x_k_smooth=NULL;
  
  long int t0=clock();
  double model_loglik=MISSING_VALUE, 
    model_dic1, eff_num_param1, model_dic2, eff_num_param2;

  params *pars=NULL;
  double *lliks=NULL, *lpriors=NULL;
  int num_lliks=0, num_lpriors=0;
  
  if(inpars_numpars>0)
    loglik(NULL); // gives us the global 'numpar', 

  // DEBUG 
  Rcout << "inpars_numpars=" << inpars_numpars << " numpar=" << numpar << 
    " inpars_numsets=" << inpars_numsets << " no_inpars=" << no_inpars <<
    std::endl;
  
  if(no_inpars)
    {
      // DEBUG Rcout << "Doing MCMC" << std::endl;
      pars=layer_mcmc(num_mcmc, burnin, spacing, numtemp, 
		      do_importance ? true : false, 
		      dosmooth | do_realizations, do_realizations,
		      NULL, &x, NULL, T_ground, 
		      &model_loglik, 
		      &model_dic1, &eff_num_param1, 
		      &model_dic2, &eff_num_param2);

      if(return_mcmc)
	{
	  num_lliks=num_lpriors=num_mcmc;
	  lliks=new double[num_mcmc];
	  lpriors=new double[num_mcmc];
	  for(j=0;(int)j<num_mcmc;j++)
	    {
	      lliks[j]=pars[j].log_lik;
	      lpriors[j]=pars[j].log_prob-pars[j].log_lik;
	    }
	}
    }
  else if(inpars_numsets>=1 && !dosmooth && !do_realizations)
    {
      Rcout << "*********************************" << std::endl;
      Rcout << "*********************************" << std::endl;
      Rcout << " Starting loglik treatment" << std::endl;
      Rcout << "*********************************" << std::endl;
      Rcout << "*********************************" << std::endl;
	
      if(!silent && inpars_numsets==1)
	detailed_loglik=true;
      
      num_lliks=num_lpriors=inpars_numsets;
      lliks=new double[inpars_numsets];
      lpriors=new double[inpars_numsets];
      
      for(int k=0;k<inpars_numsets;k++)
	{
	  int num_missing=0;
	  for(j=0;(int)j<inpars_numpars;j++)
	    if(in_pars(k,j)==MISSING_VALUE)
	      num_missing++;
	    
	  pars=new params[1];
	  pars[0].numparam=inpars_numpars;
	  pars[0].param=new double[inpars_numpars];

	  Rcout << "k=" << k << std::endl;
	
	  if(num_missing>0)
	    {
	      Rcout << "optim started" << std::endl;

	      bool detailed_loglik_buffer=detailed_loglik;
	      double lbest=-1e+200;
	      NumericVector outpars(inpars_numpars);
	      
	      detailed_loglik=false;
	      for(i=0;(int)i<num_optim;i++)
		{
		  NumericVector initpars(num_missing);
		  int jj=0;
		  for(j=0;(int)j<inpars_numpars;j++)
		    if(in_pars(k,j)==MISSING_VALUE)
		      {
			initpars[jj]=init_par(j);
			jj++;
		      }
		  partial_par=new double[inpars_numpars];
		  for(j=0;(int)j<inpars_numpars;j++)
		    partial_par[j]=transform_parameter(in_pars(k,j),
						       par_trans_type[j]);
		  Environment stats("package:stats");
		  Function optim=stats["optim"];
		  List optres=optim(Rcpp::_["par"]= initpars,
				    Rcpp::_["fn"]=Rcpp::InternalFunction(partial_loglikwrapper),
				    Rcpp::_["method"]="BFGS");
		  NumericVector outpars2=as<NumericVector>(optres["par"]);
		  NumericVector outval=as<NumericVector>(optres["value"]);
		  
		  double ll=-outval(0);
		  if(ll > lbest)
		    {
		      lbest=-outval(0);
		      jj=0;
		      for(j=0;(int)j<inpars_numpars;j++)
			{
			  if(in_pars(k,j)==MISSING_VALUE)
			    {
			      outpars(j)=outpars2(jj);
			      jj++;
			    }
			  else
			    {
			      outpars(j)=in_pars(k,j);
			    }
			}
		    }
		}
	      
	      detailed_loglik=detailed_loglik_buffer;
	      Rcout << "optim loop done" << std::endl;

	      for(j=0;(int)j<inpars_numpars;j++)
		{
		  if(in_pars(k,j)!=MISSING_VALUE)
		    {
		      pars[0].param[j]=transform_parameter(in_pars(k,j),
							   par_trans_type[j]);
		    }
		  else
		    {
		      pars[0].param[j]=outpars(j);
		    }
		}
	      
	      Rcout << "fill out done" << std::endl;
	
	      delete [] partial_par;
	    }
	  else
	    {
	      Rcout << "fill out started" << std::endl;
	      for(j=0;(int)j<inpars_numpars;j++)
		pars[0].param[j]=transform_parameter(in_pars(k,j),
						     par_trans_type[j]);
	      Rcout << "fill out done" << std::endl;
	    }

	  show_vec("params",pars[0].param,inpars_numpars);
	  Rcout << "Fill in loglik and logprior" << std::endl;
	  lliks[k]=loglik(pars[0].param, 0, 0);
	  Rcout << "Filled in loglik" << std::endl;
	  lpriors[k]=logprob(pars[0],0.0)-lliks[k];
	  Rcout << "Filled in logprior" << std::endl;
	}
      
      Rcout << " Ending loglik treatment" << std::endl;
    }
  else if(inpars_numsets==1 && dosmooth && num_smooth==1) 
    // parameter-estimate-based smoothing
    {
      // PS: Should *not* have any missing values!
      // DEBUG: Rcout << "Fetching smoothing for estimates" << std::endl;

      num_mcmc=inpars_numsets;
      pars=new params[1];
      pars[0].numparam=numpar;
      pars[0].param=new double[numpar];
      for(i=0;i<numpar;i++)
	pars[0].param[i]=in_pars(0,i);
      loglik(pars[0].param, 1, 0);

      P_k_smooth=new double **[meas_smooth_len];
      x_k_smooth=new double *[meas_smooth_len];

      for(unsigned int k=0;k<meas_smooth_len;k++)
	{
	  x_k_smooth[k]=new double[num_states];
	  P_k_smooth[k]=Make_matrix(num_states,num_states);

	  for(i=0;i<num_states;i++)
	    {
	      x_k_smooth[k][i]=x_k_s_kept[k][i];
	      for(j=0;j<(int)num_states;j++)
		P_k_smooth[k][i][j]=P_k_s_kept[k][i][j];
	    }
	}
      
      cleanup_x_and_P(meas_smooth_len);
    }
  else if(inpars_numsets>=1 && dosmooth)
    {
      // DEBUG: Rcout << "Getting MCMC from in_pars" << std::endl;
      
      num_mcmc=inpars_numsets;
      pars=new params[num_mcmc];
      for(j=0;(int)j<num_mcmc;j++)
	{
	  pars[j].numparam=numpar;
	  pars[j].param=new double[numpar];
	  for(i=0;i<numpar;i++)
	    pars[j].param[i]=in_pars(j,i);
	}
      // DEBUG: Rcout << "Done getting MCMC from in_pars" << std::endl;

      x=new double**[num_states];
      for(i=0;i<num_states;i++)
	{
	  x[i]=new double*[meas_smooth_len];
	  for(unsigned int k=0;k<meas_smooth_len;k++)
	    {
	      x[i][k]=new double[num_smooth*num_mcmc];
	      for(j=0;j<(int)(num_smooth*num_mcmc);j++)
		x[i][k][j]=0.0;
	    }
	}

      for(i=0;(int)i<num_mcmc;i++)
	{
	  loglik(pars[i].param, 1, 0);

	  for(unsigned int k=0;k<meas_smooth_len;k++)
	    {
	      double *xx=new double[num_states], 
		**ss=Make_matrix(num_states,num_states);
	      unsigned int ii,jj;
	      
	      for(ii=0;ii<num_states;ii++)
		xx[ii]=x_k_s_kept[k][ii];
	      for(ii=0;ii<num_states;ii++)
		for(jj=0;jj<num_states;jj++)
		  ss[ii][jj]=P_k_s_kept[k][ii][jj];
		  
	      int is_singular=0;
	      

	      double det=matrix_determinant(ss,num_states);
	      double diagprod=1.0;
	      for(ii=0;ii<num_states;ii++)
		diagprod*=ss[ii][ii];
	      is_singular=(det <= ABSVAL((diagprod*1e-8)) ? 1 : 0);


	      // sample:
	      double **new_x;
	      if(!is_singular || num_states==1)
		//new_x=sample_from_multinormal(num_smooth,xx,ss,
		//			      num_states,rptr);
		new_x=multinormal_sample(num_smooth,xx,ss,
					 num_states);
	      else
		{
		  new_x=new double*[num_smooth];
		  for(ii=0;ii<num_smooth;ii++)
		    {
		      new_x[ii]=new double[num_states];
		      double new_x_diag=sqrt(ss[0][0])*get_random_gauss();
		      
		      for(jj=0;jj<num_states;jj++)
			{
			  new_x[ii][jj]=xx[jj]+new_x_diag;
			}
		    }
		}

	      for(j=0;j<(int)num_smooth;j++)
		{
		  for(ii=0;ii<num_states;ii++)
		    x[ii][k][i*num_smooth+j]=new_x[j][ii];
		}

	      doubledelete(new_x,num_smooth);
	      doubledelete(ss,numsites);
	      delete [] xx;
	    }

	  cleanup_x_and_P(meas_smooth_len);
	}
    }
  long int t1=clock();

  
    
  
  unsigned int numsamples=(unsigned int) num_mcmc;
  // DEBUG:
  Rcout << "numsamples=" << numsamples << std::endl;


  // Something goes wrong here, for loglik mode!
  double **parsample=Make_matrix(numpar,numsamples);
  double **parsample_repar=Make_matrix(numpar,numsamples);
  if(no_inpars || dosmooth || do_realizations)
    {
      for(i=0;i<numpar;i++)
	for(j=0;j<(int)numsamples;j++)
	  {
	    parsample[i][j]=invtransform_parameter(pars[j].param[i], par_trans_type[i]);
	    parsample_repar[i][j]=pars[j].param[i];
	  }
    }
  
  // DEBUG:
  Rcout << " parameter samples fetched" << std::endl;

  char **par_name_orig=new char*[numpar];
  char **par_name_new=new char*[numpar];
  for(i=0;i<numpar;i++)
    {
      par_name_orig[i]=new char[300];
      strcpy(par_name_orig[i],par_name[i]);
    }
  
  // DEBUG:
  Rcout << " parameter names fetched" << std::endl;

  // transform the expectancies from logarithmic size to
  // original size:
  for(i=0;i<numpar;i++)
    {
      s=par_series[i];
      par_name_new[i]=new char[300];
      strcpy(par_name_new[i],par_name_orig[i]);
      if(ser[s].pr->is_log)
	{
	  if(par_type[i]==MU)
	    {
	      snprintf(par_name_new[i],99,"exp(%s)",par_name_orig[i]);
	      if(no_inpars || dosmooth || do_realizations)
		for(j=0;j<(int)numsamples;j++)
		  parsample[i][j]=exp(parsample[i][j]);
	    }
	  if(par_type[i]==OBS_SD)
	    {
	      snprintf(par_name_new[i],99,"%s_origscale",par_name_orig[i]);
	      if(no_inpars || dosmooth || do_realizations)
		for(j=0;j<(int)numsamples;j++)
		  parsample[i][j]=(ser[s].pr->is_log==2 ? 
				   exp(ser[s].mean_val) : ser[s].mean_val)*
		    sqrt(exp(parsample[i][j]*parsample[i][j])-1.0)*
		    exp(parsample[i][j]*parsample[i][j]/2.0);
	    }
	}
    }

  for(i=0;i<numpar;i++)
    Rcout << i << " " << par_name_orig[i] << " " << par_name_new[i] << std::endl;
  
  
  // DEBUG:
  Rcout << "Stationary stdev?" << std::endl;

  if(use_stationary_sdev)
  if(no_inpars || dosmooth || do_realizations)
    {
      for(s=0;s<num_series;s++)
	{
	  for(unsigned int l=0;l<ser[s].numlayers;l++)
	    {
	      bool all_global=false;
	      for(j=0;j<(int)numsites && !all_global;j++)
		{
		  int  index_diffusion=-1, index_pull=-1;
		  bool found_diffusion=false, found_pull=false;
		  bool global_diffusion=false, global_pull=false;
		  
		  for(i=0;i<numpar && (!found_diffusion || !found_pull);i++)
		    {
		      if(par_series[i]==(int)s && par_layer[i]==(int)(l+1) &&
			 (par_region[i]==(int)j || par_region[i]<0))
			{
			  if(par_type[i]==SIGMA)
			    {
			      found_diffusion=true;
			      index_diffusion=i;
			      if(par_region[i]<0)
				global_diffusion=true;
			    }
			  else if(par_type[i]==DT)
			    {
			      found_pull=true;
			      index_pull=i;
			      if(par_region[i]<0)
				global_pull=true;
			    }
			}
		    }
		  
		  if(global_diffusion && global_pull)
		    all_global=true;
		  
		  if(found_diffusion && found_pull)
		    {
		      char s_name[1000];
		      
		      if(all_global)
			snprintf(s_name, 999,"sd_%s_%d", ser[s].name,l+1);
		      else
			snprintf(s_name, 999,"sd_%s_%d_r%d", ser[s].name,l+1,j);
		      strcpy(par_name_new[index_diffusion],s_name);
		      
		      for(i=0;i<numsamples;i++)
			parsample[index_diffusion][i]=
			  parsample[index_diffusion][i]*
			  sqrt(parsample[index_pull][i]/2.0);
		    }
		  else
		    {
		      if(!found_pull)
			{
			  Rcout << "Stationary standard deviation for series " << s <<
			    " layer " << l+1 << " site " << j << std::endl;
			  Rcout << "could not be calculated due to no pull given "
			    "for this process " << std::endl;
			  Rcout << "(i.e. the process is not stationary)." << std::endl;
			  //else
			  //Rcout << "Couldn't find diffusion for series " << s << 
			  //  " layer " << l+1 << " site " << j << "!" << std::endl;
			}
		    }
		}
	    }
	}
    }
  
  // DEBUG:
  Rcout << "Add extra cycle info" << std::endl;

  int numpar2=numpar;
  bool has_cycles[10];
  if(num_series_feed>0)
    {
      numpar2++;
      for(i=0;i<10;i++)
	{
	  has_cycles[i]=false;
	  for(j=0;j<(int)numsamples && !has_cycles[i];j++)
	    if(pars[j].cycles[i]!=MISSING_VALUE && (ABSVAL((pars[j].cycles[i]))>1e-30))
	      {
                if(!silent)
		  Rcout << "j=" << j << " i=" << i << " cycles=" << 
		    pars[j].cycles[i] << std::endl;
		has_cycles[i]=true;
	      }
	  if(has_cycles[i])
	    numpar2++;
	}
    }
  StringVector parnames(numpar2);
  
  // DEBUG:
  Rcout << "best pars repar" << std::endl;

  double *best_pars_repar=new double[numpar];
  if(no_inpars || dosmooth || do_realizations)
    if(!do_ml)
      for(j=0;j<(int)numpar;j++)
	best_pars_repar[j]=
	  find_statistics(parsample_repar[j], numsamples, MEDIAN);
      
  // DEBUG:
  Rcout << "ML, do_ml=" << do_ml << std::endl;

  double aic=MISSING_VALUE, aicc=MISSING_VALUE, bic=MISSING_VALUE;
  double *ml_pars=new double[numpar];
  double best_loglik=MISSING_VALUE;
  if(no_inpars || dosmooth || do_realizations)
    {
  if(do_ml)
    {
      // traverse the number of wanted hill-climbs:
      
      if(!silent)
	Rcout << "ML started" << std::endl;
      ml_started=true;
      
      for(i=0;(int)i<num_optim;i++)
	{
	  if(!silent)
	    Rcout << "Optimization nr. " << i << std::endl;
	  
	  double *curr_par=new double[numpar];
	  
	  for(j=0;j<(int)numpar;j++)
	    {
	      if(i==0)
		curr_par[j]=find_statistics(parsample_repar[j],numsamples,MEDIAN);
	      else if(i==1)
		curr_par[j]=find_statistics(parsample_repar[j],numsamples,MEAN);
	      else if(num_optim==3)
		curr_par[j]=parsample_repar[j][numsamples/2];
	      else
		curr_par[j]=parsample_repar[j][(i-2)*(numsamples-1)/(num_optim-3)];
	    }
	  
	  
	  
	  /*
	  // do the optimization:
	  int maxiter=1000;
	  //double *pars2=gsl_optimization_cover(minusloglik, 
	  double *pars2=quasi_newton(minusloglik,
	  numpar, // number of parameters 
	  curr_par, // starting values  
	  0.00000001, // precision 
	  maxiter, // Maximum number of
	  // iterations. The number of iterations
	  // neccessary will be stored here after the
	  // optimization and can be checked 
	  true,false,
	  silent);
	  */
	  
	  
	  NumericVector initpars(numpar);
	  for(j=0;j<(int)numpar;j++)
	    initpars[j]=curr_par[j];
	  Environment stats("package:stats");
	  Function optim=stats["optim"];
	  List optres=optim(Rcpp::_["par"]= initpars,
			    Rcpp::_["fn"]=Rcpp::InternalFunction(loglikwrapper),
			    Rcpp::_["method"]="BFGS");
	  NumericVector outpars=as<NumericVector>(optres["par"]);
	  double *pars2=new double[numpar];
	  for(j=0;j<(int)numpar;j++)
	    pars2[j]=outpars[j];
	  
	  
	  /* RInside R;
	     R["loglikwrapper"] = Rcpp::InternalFunction(&loglikwrapper);
	     R["initpars"] = initpars;
	     
	     NumericVector initpars(numpar);
	     for(j=0;j<numpar;j++)
	     initpars[j]=curr_par[j];
	     List optres=R.parseEval("optim(initpars,loglikwrapper)");
	     NumericVector outpars=as<NumericVector>(optres["par"]);
	     double *pars2=new double[numpar];
	     for(j=0;j<numpar;j++)
	     pars2[j]=outpars[j];
	  */
	  
	  
	  // Fetch the log-likelihood for the optimized parameters:
	  double curr_loglik=loglik(pars2);
	  
	  // Check if this is the best result so far:
	  if(i==0 || curr_loglik>best_loglik)
	    {
	      // If so, store the parameter array and the log-likelihood:
	      for(j=0;j<(int)numpar;j++)
		best_pars_repar[j]=pars2[j];
	      best_loglik=curr_loglik;
	    }
	  
	  delete [] curr_par;
	  delete [] pars2;
	}
      
      // transform the expectancies from logarithmic size to
      // original size:
      for(j=0;j<(int)numpar;j++)
	{
	  ml_pars[j]=invtransform_parameter(best_pars_repar[j], 
					    par_trans_type[j]);
	      
	  s=par_series[j];
	  if(ser[s].pr->is_log)
	    {
	      if(par_type[j]==MU)
		ml_pars[j]=exp(ml_pars[j]);
	      if(par_type[j]==OBS_SD)
		ml_pars[j]=ser[s].mean_val*
		  sqrt(exp(ml_pars[j]*ml_pars[j])-1.0)*
		  exp(ml_pars[j]*ml_pars[j]/2.0);
	    }
	}
      
      if(return_residuals==1)
	{
	  if(!silent)
	    Rcout << "Running loglik to fetch ML residuals..." << std::endl;
	  
	  double lres=loglik(best_pars_repar, 0, 0, 0, NULL, 0, NULL,
			     1, &resids_time,
				 &resids, &prior_expected_values,
			     &resid_numcol, &resid_len);
	  
	  if(!resids)
	    Rcout << "Running loglik to fetch ML residuals failed! lres=" <<
		  lres << std::endl;
	  else if(!silent)
	    Rcout << "lres=" << lres << std::endl;
	}
      //loglik(best_pars_repar, 0, 0,1,filestart);
      
      
      
      int k=numpar;
      int nn=0;
      for(s=0;s<num_series;s++)
	nn+=double(ser[s].meas_len);
      
      aic=-2*best_loglik + 2.0*double(k);
      aicc=-2*best_loglik + 2.0*double(k) + double(k*(k+1))/double(nn-k-1);
      bic=-2*best_loglik + log(double(nn))*double(k);
    }
  if(!do_ml && return_residuals==1)
    {
      double *medpar=new double[numpar];
      for(j=0;j<(int)numpar;j++)
	medpar[j]=find_statistics(parsample_repar[j], numsamples, MEDIAN);
      
      if(!silent)
	Rcout << "Running loglik to fetch Bayesian median residuals..." << std::endl;
	  
      double lres=loglik(medpar, 0, 0, 0, NULL, 0, NULL,
			 1, &resids_time, &resids, &prior_expected_values,
			 &resid_numcol, &resid_len);
      
      if(!resids)
	Rcout << "Running loglik to fetch ML residuals failed! lres=" <<
	  lres << std::endl;
      else if(!silent)
	Rcout << "lres=" << lres << std::endl;
	  
      delete [] medpar;
    }
    }
  
  Rcout << "ML done" << std::endl;

#ifdef DETAILED_TIMERS
  List alltimers=List::create(
		     Named("loglik.time",double(timers[1][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("loglik.init.time",double(timers[2][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("loglik.mainloop.time",double(timers[3][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("loglik.smoother.time",double(timers[4][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("loglik.realization.time",double(timers[5][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("logprior.time",double(timers[9][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("matrix.inverse.time",double(timers[10][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("complex.matrix.inverse.time",double(timers[11][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("matrix.eigenvectors.time",double(timers[12][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("multinormal.sample.time",double(timers[13][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("pdf.multinormal.time",double(timers[14][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("variance.estimation.time",double(timers[15][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("mcmc.init.time",double(timers[20][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("mcmc.burnin.time",double(timers[21][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("mcmc.mainloop.time",double(timers[22][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("model.fit.estimation.time",double(timers[23][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("test.ou.time",double(timers[24][2])/
			   double(CLOCKS_PER_SEC)),
		     Named("runs.test.time",double(timers[25][2])/
		     	   double(CLOCKS_PER_SEC)),
		     Named("analyze.residuals.time",double(timers[26][2])/
			   double(CLOCKS_PER_SEC)));
#endif // DETAILED_TIMERS
  
  
  List out;
  if(!do_ml)
    out=List::create(Named("model.log.lik",model_loglik),
		     Named("computing.time", double(t1-t0)/
			   double(CLOCKS_PER_SEC)),
#ifdef DETAILED_TIMERS
		     Named("all.timers",alltimers),
#endif // DETAILED_TIMERS
		     Named("model.dic1",model_dic1),
		     Named("eff.num.param1",eff_num_param1),
		     Named("model.dic2",model_dic2),
		     Named("eff.num.param2",eff_num_param2));
  else
    out=List::create(Named("model.log.lik",model_loglik),
		     Named("computing.time", double(t1-t0)/
			   double(CLOCKS_PER_SEC)),
#ifdef DETAILED_TIMERS
		     Named("all.timers",alltimers),
#endif // DETAILED_TIMERS
		     Named("model.dic1",model_dic1),
		     Named("eff.num.param1",eff_num_param1),
		     Named("model.dic2",model_dic2),
		     Named("eff.num.param2",eff_num_param2),
		     Named("ML.loglik",best_loglik),
		     Named("AIC",aic),
		     Named("AICc",aicc),
		     Named("BIC",bic));
  
  // DEBUG: Rcout << "Start summary statistics (1)" << std::endl;

  // insert mean, median and 95% credibility interval for
  // each parameter into the return list:
  if(!silent)
    Rcout << "numpar=" << numpar << std::endl;

  double *est_pars=new double[numpar];
  if(no_inpars || dosmooth || do_realizations)
    {
      for(j=0;j<(int)numpar;j++)
	{
	  double mean= find_statistics(parsample[j], numsamples, MEAN);
	  double med=  find_statistics(parsample[j], numsamples, MEDIAN);
	  double lower=find_statistics(parsample[j], numsamples, PERCENTILE_2_5);
	  double upper=find_statistics(parsample[j], numsamples, PERCENTILE_97_5);
	  
	  if(use_half_times && !strncmp(par_name[j],"dt_",3))
	    {
	      mean*=log(2.0);
	      med*=log(2.0);
	      lower*=log(2.0);
	      upper*=log(2.0);
	      if(do_ml && ml_pars)
		ml_pars[j]*=log(2.0);
	    }
	  
	  List parlist;
	  if(!do_ml)
	    {
	      est_pars[j]=mean;
	      parlist=List::create(Named("mean",mean),
				   Named("median", med),
				   Named("lower02_5",lower),
				   Named("upper97_5",upper) );
	    }
	  else
	    {
	      est_pars[j]=ml_pars[j];
	      parlist=List::create(Named("ML",ml_pars[j]),
				   Named("mean",mean),
				   Named("median", med),
				   Named("lower02_5",lower),
				   Named("upper97_5",upper) );
	    }
	  out[par_name_new[j]]=parlist;
	  
	  if(!silent)
	    {
	      Rcout << "j=" << j << std::endl;
	      Rcout << par_name_new[j] << std::endl;
	      Rcout << "mean=" << mean << " median=" << med << 
	    " lower=" << lower << " upper=" << upper << std::endl;
	    }
	  
	  std::string pname(par_name_new[j]);
	  parnames(j)=pname;
	}
    }
  if(ml_pars)
    delete [] ml_pars;

  if(no_inpars || dosmooth || do_realizations)
    {
      NumericVector est_paramset(numpar);
      for(i=0;i<numpar;i++)
	est_paramset(i)=est_pars[i];
      out["est.par"]=est_paramset;
      delete [] est_pars;
    }
  

  
  // DEBUG:
  Rcout << "Start summary statistics (2)" << std::endl;

  if(no_inpars || dosmooth || do_realizations)
  if(num_series_feed>0)
    {
      double *iscomplex=new double[numsamples];
      int k=numpar;
      
      for(j=0;j<(int)numsamples;j++)
	iscomplex[j]=pars[j].iscomplex ? 1.0 : 0.0;
      
      double ic_mean=find_statistics(iscomplex, numsamples, MEAN);
      double ic_med=find_statistics(iscomplex, numsamples, MEDIAN);
      double ic_lower=find_statistics(iscomplex, numsamples, PERCENTILE_2_5);
      double ic_upper=find_statistics(iscomplex, numsamples, PERCENTILE_97_5);
      
      delete [] iscomplex;
      
      List parlist=List::create(Named("mean",ic_mean),
				Named("median", ic_med),
				Named("lower02_5",ic_lower),
				Named("upper97_5",ic_upper) );
      
      out["complex.eigen"]=parlist;
      
      std::string ic_name("complex.eigen");
      parnames(k)=ic_name;
      k++;
      
      for(i=0;i<10;i++)
	if(has_cycles[i])
	  {
	    double *cycle=new double[numsamples];
	    for(j=0;j<(int)numsamples;j++)
	      cycle[j]=pars[j].cycles[i];
	    
	    double c_mean=find_statistics(cycle, numsamples, MEAN,0.01);
	    double c_med=find_statistics(cycle, numsamples, MEDIAN,0.01);
	    double c_lower=find_statistics(cycle, numsamples, PERCENTILE_2_5,0.01);
	    double c_upper=find_statistics(cycle, numsamples, PERCENTILE_97_5,0.01);
	    
	    delete [] cycle;
	    
	    List parlist=List::create(Named("mean",c_mean),
				      Named("median", c_med),
				      Named("lower02_5",c_lower),
				      Named("upper97_5",c_upper) );
	    char cyclename[100];
	    snprintf(cyclename,99,"cycle%02d",i+1);
	    out[cyclename]=parlist;
	    
	    std::string c_name(cyclename);
	    parnames(k)=c_name;
	    k++;
	  }
    }
  
  // DEBUG:
  Rcout << "Set parnames" << std::endl;
  
  out["parameter.names"]=parnames;

  // DEBUG:
  Rcout << "Residuals" << std::endl;

  if(no_inpars || dosmooth || do_realizations)
  if(return_residuals==1)
    {
      if(resids && resid_len>0)
	{
	  NumericMatrix standardized_residuals(resid_len, resid_numcol);
	  for(i=0;(int)i<resid_numcol;i++)
	    for(j=0;(int)j<resid_len;j++)
	      standardized_residuals(j,i)=resids[i][j];
	  out["standardized.residuals"]=standardized_residuals;
	  doubledelete(resids, resid_numcol);
	}
      else
	{
	  Rcout << "No residuals found!" << std::endl;
	}
      
      if(resids_time && resid_len>0)
	{
	  NumericVector residuals_time(resid_len);
	  for(j=0;(int)j<resid_len;j++)
	    residuals_time(j)=resids_time[j];
	  out["residuals.time"]=residuals_time;
	  delete [] resids_time;
	}

      if(prior_expected_values && resid_len>0)
	{
	  NumericMatrix prior_expectancy(resid_len, resid_numcol);
	  for(i=0;(int)i<resid_numcol;i++)
	    for(j=0;(int)j<resid_len;j++)
	      prior_expectancy(j,i)=prior_expected_values[i][j];
	  out["prior.expected.values"]=prior_expectancy;
	  doubledelete(prior_expected_values, resid_numcol);
	}
    }
   
  // DEBUG:
  Rcout << "Start smoothing" << std::endl;

  if(dosmooth && !do_realizations)
    {
      StringVector process_names(num_states);
      NumericVector process_time_points(meas_smooth_len);
      NumericMatrix process_mean(num_states,meas_smooth_len);
      NumericMatrix process_lower95(num_states,meas_smooth_len);
      NumericMatrix process_upper95(num_states,meas_smooth_len);
      char **procname=new char*[num_states];
      unsigned int k;
      
      for(j=0;j<(int)meas_smooth_len;j++)
	process_time_points(j)=meas_smooth[j].tm;
      
      for(i=0;i<num_states;i++)
	{
	  int s=state_series[i], l=state_layer[i];
	  int site=i%numsites;
	  
	  procname[i]=new char[1000];
	  if(numsites<=1)
	    snprintf(procname[i], 999,"%s_layer%03d",ser[s].name,l+1);
	  else
	    snprintf(procname[i], 999,
		     "%s_site%03d_layer%03d",ser[s].name,site,l+1);
	  
	  std::string pname(procname[i]);
	  
          process_names(i)=pname;
	  
	  if(inpars_numsets==1 && num_smooth==1) 
	    // parameter-estimate-based smoothing
	    {
	      // DEBUG: Rcout << "Making standard smoothing output" << std::endl;
	      for(k=0;k<meas_smooth_len;k++)
		{
		  process_mean(i,k)=x_k_smooth[k][i];
		  process_lower95(i,k)=x_k_smooth[k][i]-
		    1.96*sqrt(P_k_smooth[k][i][i]);
		  process_upper95(i,k)=x_k_smooth[k][i]+
		    1.96*sqrt(P_k_smooth[k][i][i]);
		}
	    }
	  else
	    {
	      for(k=0;k<meas_smooth_len;k++)
		{
		  process_mean(i,k)=find_statistics(x[i][k], num_smooth*numsamples, MEAN);
		  process_lower95(i,k)=
		    find_statistics(x[i][k], num_smooth*numsamples, PERCENTILE_2_5);
		  process_upper95(i,k)=
		    find_statistics(x[i][k], num_smooth*numsamples, PERCENTILE_97_5);
		}
	    }
	}
      
      out["process.names"]=process_names;
      out["process.time.points"]=process_time_points;
      out["process.mean"]=process_mean;
      out["process.lower95"]=process_lower95;
      out["process.upper95"]=process_upper95;
      
     
      if(do_return_smoothing_samples)
	{
	  NumericMatrix firstmatrix(meas_smooth_len,num_smooth*numsamples);
	  for(j=0;j<(int)meas_smooth_len;j++)
	    for(k=0;k<num_smooth*numsamples;k++)
	      firstmatrix(j,k)=x[0][j][k];
	  List smooth_samples=List::create(Named(procname[0],firstmatrix));
	  
	  for(i=1;i<num_states;i++)
	    {
	      NumericMatrix nextmatrix(meas_smooth_len,num_smooth*numsamples);
	      for(j=0;j<(int)meas_smooth_len;j++)
		for(k=0;k<num_smooth*numsamples;k++)
		  nextmatrix(j,k)=x[i][j][k];
	      smooth_samples[procname[i]]=nextmatrix;
	    }
	  
	  out["smoothing.samples"]=smooth_samples;
	}
      
      if(inpars_numsets==1 && num_smooth==1) 
	// parameter-estimate-based smoothing
	{
	  // DEBUG: Rcout << "Making P_k smoothing output" << std::endl;
	  // Rcout << "P_k_smooth" << (P_k_smooth ? 1 : 0) << std::endl;
	  NumericMatrix firstmatrix(num_states, num_states);
	  for(i=0;i<num_states;i++)
	    for(j=0;j<(int)num_states;j++)
	      firstmatrix(i,j)=P_k_smooth[0][i][j];
	  char strbuff[20]="P_k_00001";
	  List P_k=List::create(Named(strbuff,firstmatrix));
	  
	  for(k=1;k<meas_smooth_len;k++)
	    {
	      NumericMatrix nextmatrix(num_states, num_states);
	      for(i=0;i<num_states;i++)
		for(j=0;j<(int)num_states;j++)
		  nextmatrix(i,j)=P_k_smooth[k][i][j];
	      snprintf(strbuff, 19, "P_k_%05d",k+1);
	      P_k[strbuff]=nextmatrix;
	    }
	  
	  out["P.k.smooth"]=P_k;
	}

      doubledelete(procname,num_states);
    }
  
  if(x)
    tripledelete(x,num_states,meas_smooth_len);
  
  if(x_k_smooth)
    doubledelete(x_k_smooth,meas_smooth_len);
  
  if(P_k_smooth)
    tripledelete(P_k_smooth,meas_smooth_len,num_states);
  
  
  // DEBUG:
  Rcout << "Smoothing done" << std::endl;

  if(do_realizations && numit_realizations_made>0 && x_k_realized_all!=NULL)
    {
      NumericVector realization_time_points(meas_smooth_len);
      for(j=0;j<(int)meas_smooth_len;j++)
	realization_time_points(j)=meas_smooth[j].tm;
      out["realization.time.points"]=realization_time_points;
      
      char **procname=new char*[num_states];
      for(i=0;i<num_states;i++)
	{
	  int s=state_series[i], l=state_layer[i];
	  int site=i%numsites;
	  
	  procname[i]=new char[200];
	  if(numsites<=1)
	    snprintf(procname[i], 199,"%s_layer%03d",ser[s].name,l+1);
	  else
	    snprintf(procname[i], 199,
		     "%s_site%03d_layer%03d",ser[s].name,site,l+1);
	}
      
      unsigned int k;
      NumericMatrix firstmatrix(meas_smooth_len,numit_realizations_made);
      for(j=0;j<(int)meas_smooth_len;j++)
	for(k=0;k<numit_realizations_made;k++)
	  {
	    if(!silent)
	      Rcout << "j=" << j << " k=" << k << " x[k]=" << x_k_realized_all[k] << " x_k_j=" << x_k_realized_all[k][j] << std::endl;
	    firstmatrix(j,k)=x_k_realized_all[k][j][0];
	  }
      List realizations=List::create(Named(procname[0],firstmatrix));
      
      for(i=1;i<num_states;i++)
	{
	  // DEBUG: Rcout << "state=" << i << std::endl;
	  NumericMatrix nextmatrix(meas_smooth_len,numit_realizations_made);
	  for(j=0;j<(int)meas_smooth_len;j++)
	    for(k=0;k<numit_realizations_made;k++)
	      {
		if(!silent)
		  Rcout << "state=" << i << " j=" << j << " k=" << k << " x[k]=" << x_k_realized_all[k] << " x_k_j=" << x_k_realized_all[k][j] << std::endl;
		nextmatrix(j,k)=x_k_realized_all[k][j][i];
	      }
	  // DEBUG: Rcout << "Set realizations number " << i << std::endl;
	  realizations[procname[i]]=nextmatrix;
	}
      
      if(!silent)
	Rcout << "Set realizations" << std::endl;
      out["realizations"]=realizations;
      
      if(!silent)
	Rcout << "Cleanup" << std::endl;
      doubledelete(procname,num_states);
    }
  
  // DEBUG:
  Rcout << "Realizations done" << std::endl;
  
  if(no_inpars || dosmooth || do_realizations)
    {
      NumericVector est_origpar(numpar);
      for(i=0;i<numpar;i++)
	est_origpar(i)=best_pars_repar[i];
      out["est.origpar"]=est_origpar;
      delete [] best_pars_repar;
    }
  
  // DEBUG:
  Rcout << " Returning MCMC" << std::endl;

  if(return_mcmc)
    {
      NumericMatrix mcmcsamples(numpar,numsamples);
      for(i=0;i<numpar;i++)
	for(j=0;j<(int)numsamples;j++)
	  mcmcsamples(i,j)=parsample[i][j];
      out["mcmc"]=mcmcsamples;
      
      NumericMatrix mcmcsamples2(numpar,numsamples);
      for(i=0;i<numpar;i++)
	for(j=0;j<(int)numsamples;j++)
	  mcmcsamples2(i,j)=parsample_repar[i][j];
      out["mcmc.origpar"]=mcmcsamples2;
    }

  if(lliks)
    {
      NumericVector llik(num_lliks);
      for(j=0;(int)j<num_lliks;j++)
	llik(j)=lliks[j];
      out["loglik"]=llik;
      delete [] lliks;
    }
  
  if(lpriors && !do_ml)
    {
      NumericVector lprior(num_lpriors);
      for(j=0;(int)j<num_lpriors;j++)
	lprior(j)=lpriors[j];
      out["logprior"]=lprior;
      delete [] lpriors;
    }
  
  
  // DEBUG:
  Rcout << " Returned MCMCg" << std::endl;

  doubledelete(par_name_orig,numpar);
  Rcout << " deleted par_name_orig" << std::endl;
  doubledelete(par_name_new,numpar);
  Rcout << " deleted par_name" << std::endl;
  doubledelete(parsample,numpar);
  Rcout << " deleted parsample" << std::endl;
  doubledelete(parsample_repar,numpar);
  Rcout << " deleted parsample_repar" << std::endl;
  
  delete [] pars;
  Rcout << " deleted pars" << std::endl;
  delete [] meas_buffer;
  Rcout << " deleted meas_buffer" << std::endl;
  reset_global_variables();
  
  // DEBUG:
  Rcout << " Returning" << std::endl;

  PutRNGstate();
  return(out);
}


#endif // ifndef MAIN






#ifdef MAIN // Only for when compiling as standalone program



// ********************************************************
// Usage: Print usage and exit
// ********************************************************

void usage(void)
{
  printf("Usage: layer_analyzer [run options] <series specificiations> "
	 "<numsamples>\n"
	 "<burnin> <indep> <numtemp> \n");
  printf("or\n");
  printf("layer_analyzer [run options] -S <simulation file start> "
	 "<number of simulations> <parameter list>\n");
  printf("\n");
  printf("Series specfication syntax:\n");
  printf("  -i [series options] <series name> <series file> <prior file>\n");
  printf("Series options:\n");
  printf("         -l <num of layers> : Specifies the number of layers used\n"
	 "               (default three for the main series and one for the "
	 "supplementary series)\n");
  printf("               Should be specified before the other series options\n");
  printf("         -t: linear time dependency\n");
  printf("         -T <layer>: time integration of the layer below on the \n"
	 "               specified layer (can't be the lowest layer)\n"
	 "               PS: Automatically switches on identification\n"
	 "               prior strategy 0 (no identification restriction).\n");
  printf("         -ru : Regional mu\n");
  printf("         -rt : Regional linear time dependency\n");
  printf("         -rp <layer> : Regional pull for a given layer\n");
  printf("         -rs <layer> : Regional diffusion for a given layer\n");
  printf("         -C <layer> : Correlated diffusion for a given layer\n");
  printf("         -rC <layer> : Correlated diffusion with pairwise "
	 "correlation coefficient\n");
  printf("               PS: Inefficient numerics means that #sites<=6\n");
  printf("               PSS: Fairly unparsemoneous when #sites>=4\n");
  printf("         -np : No pull in the lowest layer layer (random walk)\n");
  printf("         -ns <layer> : No diffusion on a given layer "
	 "(should not be the bottom layer)\n");
  printf("         -1s <layer> : 1D diffusion (one dimensional noise)\n");
  printf("         -SS <layer> : indicator flag gourping of correlation "
	 "at a given layer\n");
  printf("         -SC <layer> : indicator flag removal of correlation "
	 "at a given layer\n");
  printf("         -Ss <layer> : indicator flag differentiation of diffusion\n");
  printf("         -Sp <layer> : indicator flag differentiation of pull\n");
  printf("         -Su : indicator flag differentiation of mu\n");
  printf("         -St : indicator flag differentiation of linear time dependency\n");
  printf("         -i0 : initial state treatment, with initial time equal \n");
  printf("               to the first measurement\n");
  printf("               (The default is to use the stationary distribuition as\n");
  printf("               the initial insight into the process).\n");
  printf("         -it <time>: initial state treatment, with initial time specified\n");
  printf("         -iT <datetime>: initial state treatment, with initial time \n");
  printf("               specified using datetime specification.\n");
  printf("         -is : if initial state treatment, the initial state \n");
  printf("               is assumed to be the same for all sites.\n");
  printf("         -il : if initial state treatment, the initial state \n");
  printf("               is assumed to be the same for all layers in a site.\n");
  printf("         -It <time> <init> : Initial value treatment\n");
  printf("               with the same specified initial value for all sites.\n");
  printf("         -IT <datetime> <init> : Initial value treatment\n");
  printf("               with the same specified initial value for all sites.\n");
  printf("         -U : Allows for positive (unstable) pull coefficients.\n");
  printf("               Only alloweable with initial value treatment.\n");
  printf("               The prior 95%% credibility interval for the untransformed\n");
  printf("               pulls will be set to ±1/lower boundry for\n");
  printf("               the characteristic time, using the normal distribution.\n");
  printf("               Should only be used when the pull identificaiton restriciton "
	 "is switched off.\n");
  printf("          -P <period>: Substracts a trigonometric from the data with the\n"
	 "               given period but with unknown constants before cos and sin\n"
	 "               Multiple such periodic functions can be used\n"); 
  printf("\n");
  printf("The prior file should have the format "
	 "is_log;mu1;mu2;dt1;dt2;s1;s2;lin1;lin2;beta1;beta2\n");
  printf("  indicating prior lower and upper 95%% credibility limits for all "
	 "parameter types.\n");
  printf("  is_log=0 means the data should not be log-transformed\n");
  printf("  is_log=1 means the data should be log-transformed\n");
  printf("  is_log=2 means the data already has been log-transformed\n");
  printf("  Parameter priors must relate to the log-transform series if is_log!=0\n");
  printf("\n");
  printf("Run options:\n");
  printf("         -S <simulation file start> <number of realizations> <parameter list>\n"
	 "            Simulates a set of realizations conditioned on the\n"
	 "            parameter set, times, standard deviations (if applicable),\n"
	 "            and number of observations (if applicable)\n"
	 "            *but not the observations themselves*. Used for making\n"
	 "            simulation datasets\n");
  printf("         -n : null data - ignores the data file and assumes\n"
	 "              no data. For debug purposes\n");
  printf("         -e <external serie file> : Include external time serie\n");
  printf("         -s  : Silent modus. Only shows marginal data density "
	 "and parameter estimates\n");
  printf("         -H <identification prior strategy>: Specifies how to \n"
	 "              handle the identification restrictions in the prior.\n");
  printf("              Options:\n"
	 "               0 - No identification treatment. (Default) PS: Can and\n"
	 "                   even should yield multimodal characteristic times\n"
	 "               1 - Keep upper characteristic time prior. Add\n"
	 "                   lognormally to beneath-lying characteristic times.\n"
	 "               2 - Keep lower characteristic time prior. Substract\n"
	 "                   lognormally to above-lying characteristic times.\n"
	 "               3 - Keep lower characteristic time prior. Cut\n"
	 "                   depending on that on the above-lying "
	 "characteristic times. (Recommended option)\n"
	 "               4 - Keep upper characteristic time prior. Cut\n"
	 "                   depending on that on the below-lying "
	 "characteristic times.\n");
  printf("         -o <file start> <numit> <N/A/D>: Sends process realizations "
	 "(Kalman smoother\n"
	 "            results) to a set of output files starting with the\n"
	 "            specified string, instead of showing on the screen.\n"
	 "            <numit> determines the number of iterations before \n"
	 "            giving up trying to make a realization.\n"
	 "            \"N\" gives no censoring on the realizations\n"
	 "            \"A\" censors realizations that anywhere ascends beyond the\n"
	 "            value of the previous measurement point."
	 "            \"D\" censors realizations that anywhere descends beyond the\n"
	 "            value of the previous measurement point.\n"
	 "            PS:  The last two options will affect the inference.\n"
	 "            PSS: These options are only implemented for single "
	 "data file analysis.\n");
  printf("         -b: skip importance sampling estimation of BML\n");
  printf("         -R <site> : remove site\n");
  printf("         -O <site> : only use this site\n");
  printf("         -B <age> : cut data before <age>\n");
  printf("         -A <age> : cut data after <age>\n");
  printf("         -p : pre show data\n");
  printf("         -Pr : Plot re-parametrized parameters\n");
  printf("         -Ps : Plot stationary standard deviations\n");
  printf("         -d <time diff>: Toogle smoothing analysis with the given "
	 "time resolution\n");
  printf("         -D <time diff> <start time> <end time>: Toogle smoothing analysis\n");
  printf("            with the given time resolution plus start and end time.\n");
  printf("         -G <time diff> <start time> <end time>: Toogle smoothing analysis\n");
  printf("            with the given time resolution plus start and end date-time.\n");
  printf("         -F <filestart>: File start for smoothing output and parameter "
	 "plots\n  "
	 "            (toggles sending the smoothing output to files)\n");
  printf("         -ML <number of optimizations> : "
	 "Performs ML analysis (after the\n"
	 "            Bayesian analysis) and residual dependency analysis\n");
  printf("         -E : toggles showing population standard deviation vs \n"
	 "            selection pressure (only appropriate for "
	 "evoluationary data,\n"
	 "            where the first layer is phenotype and the "
	 "second is the optimum)\n"
	 "            The difference between the inference on the "
	 "first and second\n"
	 "            layer is interpreted as selection pressure.\n");
  printf("         -C <series1> <layer> <series 2> <layer>\n");
  printf("               Injects a instantaneous correlation between one layer in\n");
  printf("               series 1 and one layer in series 2\n");
  printf("               Warning: Injecting more than one instantaneous series "
	 "correlation\n");
  printf("         -rC <series1> <layer> <series 2> <layer>\n");
  printf("               Same as option -C but with site-specific "
	 "correlation parameters\n");
  printf("               may lead to incorrectly calculated BMLs\n");
  printf("         -f <series1> <layer 1> <series 2> <layer 2>\n");
  printf("               Causal link from a specific series and layer to another "
	 "specific series and layer\n");
  printf("               that doesn't go by the layering system.\n");
  printf("               Treated with a regression parameter.\n");
  printf("         -g <series1> <layer 1> <series 2> <layer 2>\n");
  printf("               Symmetric causal links from between two specific series+layer\n");
  printf("         -c <correlation file>: Specifies the series observational\n");
  printf("               correlations for each combination of site and time. \n");
  printf("               Each line should contain:\n");
  printf("                site (unless number of sites is one), time, corr_1_2,\n"); 
  printf("                corr_1_3,...,corr_1_S,corr_2_3,...corr_2_S,...,corr_S-1,S\n");
  printf("                (a total of S*(S-1)/2 correlations) where S is the number\n");
  printf("                of series. Time is specified by floating point numbers if\n");
  printf("                that is used in the data files and date/time-format if\n");
  printf("                that is used in the data files. Do not use this option\n");
  printf("                until all input series have been specified!\n");
  printf("               If some series do not have measurements for some site+time\n");
  printf("               combinations then just put 0 there. If all but one serie\n");
  printf("               lacks measurements for a site+time, that site+time can be\n");
  printf("               omitted from the correlation file\n");
  printf("         -T <T_ground>: Set ground temperature for tempering chains.\n");
  printf("               Temperatures=T_ground^(i-1) for the i'th chain. "
	 "Default ground temperature is 2.\n");
  printf("         -tl : talkative likelihood\n");
  printf("         -tb : talkative burnin\n");
  printf("         -a <parameter list>: Set start parameters\n");
  printf("         -h <number of smoothings per iteration>. Default 10\n");
  printf("               Not meaningful without '-d', '-D' or '-G'.\n");
  printf("         -N: Toggles no sites. Reads the first column in each file as time\n");
  printf("             no matter what. Must be invoked before the options for "
	 "reading.\n");
  printf("         -2: Report half-times rather than characteristic times. Half-time="
	 "log(2)*characteristic time.\n");
  printf("         -q: Report stationary standard deivations.\n");
  exit(0);
}



void test_analysis(int num_layers_real, int num_layers_checked) 
// Assumes the file "test_<k>layer.txt" 
{
  reset_global_variables();
  ser=new series[100];
  
  double loglik=MISSING_VALUE;
  char filename[1000];
  snprintf(filename,999,"test_%dlayer.txt",num_layers_real);
  
  FILE *f=fopen(filename,"r");
  if(!f || feof(f))
    {
      printf("Couldn't find \"%s\"!\n",filename);
      exit(0);
    }
  
  double **X2=Make_matrix(1000,2);
  int n=0,m=2;
  size_t len=1000;
  char *line=new char[1000];
  getline(&line,&len,f);
  do
    {
      double x1,y1;
      char *ptr=getnextdouble(line,&x1);
      if(ptr)
	ptr=getnextdouble(ptr,&y1);
      if(ptr)
	{
	  X2[n][0]=x1;
	  X2[n][1]=y1;
	  n++;
	}
      len=1000;
      getline(&line,&len,f);
    } while(!feof(f));
  fclose(f);
  
  int num_layers=num_layers_checked;
  printf("num_layers_checked=%d\n",num_layers);
  
  unsigned int s,i,j;
  if(m!=2)
    return;
  
  num_series=1;
  ser[0].set_series(X2,n,num_layers);
  
  num_states=0;
  num_tot_layers=num_layers;
  numsites=1;
  for(int l=0;l<num_layers;l++)
    {
      state_series[num_states+l]=0;
      state_layer[num_states+l]=l;
    }
  series_state_start[0]=0;
  num_states = num_layers;
  num_tot_layers = num_layers;
  
  meas_tot_len=0;
  for(s=0;s<num_series;s++)
    {
      meas_tot_len += ser[s].meas_len;
      for(i=0;i<ser[s].meas_len;i++)
	ser[s].meas[i].index=0;
    }
  series_measurements *meas_buffer=new series_measurements[meas_tot_len];	  
  j=0;
  for(s=0;s<num_series;s++)
    for(i=0;i<ser[s].meas_len;i++)
      meas_buffer[j++]=ser[s].meas[i];
  qsort(meas_buffer,(size_t) meas_tot_len,sizeof(series_measurements), 
	compare_meas);

  meas_tot=new measurement_cluster[meas_tot_len];
  j=-1;
  for(i=0;i<meas_tot_len;i++)
    {
      if(j<0 || meas_tot[j].tm!=meas_buffer[i].tm)
	j++;
      if(meas_buffer[i].meanval!=MISSING_VALUE)
	meas_tot[j].add_measurement(meas_buffer[i]);
      else if(meas_tot[j].num_measurements==0)
	{
	  meas_tot[j].tm=meas_buffer[i].tm;
	  meas_tot[j].dt=meas_buffer[i].dt;
	}
    }
  meas_tot_len=j+1;
  delete [] meas_buffer;
  delete [] line;

  double model_loglik=MISSING_VALUE;
  params *pars=layer_mcmc(100, 1000, 10, 1, 
			  true, 0, 0, 
			  NULL, NULL, NULL, 1.15, 
			  &model_loglik);
  
  delete [] pars;
  delete [] X2;

  loglik=model_loglik;
  printf("loglik=%f\n",loglik);
}


// *******************************************************
// main-loop: Fetch user options, fetch the data
// fetch the MCMC sampling and marginal data density
// and show the results:
// *******************************************************
int main(int argc, char **argv)
{
  reset_global_variables();
  ser=new series[100];
  
  randify(); // random seed
  
  if(argc==3)
    {
      id_strategy=ID_SUB_LOWER;
      
      int numlayers_real=atoi(argv[1]);
      int numlayers_upper_checked=atoi(argv[2]);
      
      silent=1;
      for(int l=1;l<=numlayers_upper_checked;l++)
	test_analysis(numlayers_real,l);
      
      reset_global_variables();
      exit(0);
    }
  
  // indicator for whether a site is to be used or not
  // default on.
  bool *use_site=new bool[1000], plot_repar=false, 
    plot_sd=false, do_ml=false, report_sd=false, report_halflife=false;
  unsigned int i, k,l,m,s; // index variables
  int j;
  // Start and end time values, default is to include all:
  double start=MISSING_VALUE,end=MISSING_VALUE;
  int dosmooth=0, num_optim=0, sel_pres=0, preshow=0;
  double diff=0.0;
  char filestart[300]="";
  char realization_file_start[1000];
  bool do_realizations=false, do_simulations=false, 
    do_importance=true, do_dt_start_end=false;
  int numsim=0;
  double *sim_param=NULL;
  char simulation_file_start[300];
  double t_start=MISSING_VALUE, t_end=MISSING_VALUE;
  HydDateTime dt_start=NoHydDateTime, dt_end=NoHydDateTime;
  double *startpar=NULL;
  double T_ground=2.0;

  // ste the initial site usage to 'on':
  for(i=0;i<1000;i++)
    use_site[i]=true;

  // traverse the user options:
  while(argc>1 && argv[1][0]=='-')
    {
      // condition on the second letter:
      switch(argv[1][1])
	{
	case 'q':
	  report_sd=true;
	  break;
	case '2':
	  report_halflife=true;
	  break;
	case 'N':
	  nosites=true;
	  break;
	case 'T':
	  T_ground=atof(argv[2]);
	  argc--;
	  argv++;
	  break;
	case 'h':
	  num_smooth=atoi(argv[2]);
	  argc--;
	  argv++;
	  break;
	case 'a':
	  loglik(NULL);
	  startpar=new double[numpar];
	  for(i=0;i<numpar;i++)
	    {
	      startpar[i]=transform_parameter(atof(argv[2]),par_trans_type[i]);
	      argc--;
	      argv++;
	    }
	  break;
	case 'S':
	  if(num_series<=0)
	    {
	      cerr << "Simulation specification only possible if model and data\n"
		"has already been specified!" << endl;
	      exit(0);
	    }
	  loglik(NULL);
	  do_simulations=true;
	  strcpy(simulation_file_start, argv[2]);
	  argc--;
	  argv++;
	  numsim=atoi(argv[2]);
	  argc--;
	  argv++;
	  sim_param=new double[numpar];
	  for(i=0;i<numpar;i++)
	    {
	      sim_param[i]=atof(argv[2]);
	      argc--;
	      argv++;
	    }
	  break;
	case 'n':
	  nodata=true;
	  break;
	case 'H':
	  id_strategy=(INDENTIFICATION_PRIOR_HANDLING) atoi(argv[2]);
	  argc--;
	  argv++;
	  break;
	case 'o':
	  do_realizations=true;
	  strcpy(realization_file_start,argv[2]);
	  numit_realization=atoi(argv[3]);
	  if(argv[4][0]=='N')
	    real_strat=NO_CENSORING;
	  else if(argv[4][0]=='A')
	    real_strat=ASCENDING_CENSORING;
	  else if(argv[4][0]=='D')
	    real_strat=DESCENDING_CENSORING;
	  else
	    {
	      printf("Unknown option for the realization censoring strategy!\n\n");
	      usage();
	      exit(0);
	    }
	  argc-=3;
	  argv+=3;
	  break;
	case 'p':
	  preshow=1;
	  break;
	case 'i':
	  argc--;
	  argv++;
	  ser[num_series].read_series(argv,argc,num_series==0 ? true : false,
				      num_series,use_site, start,end);
	  num_series++;
	  break;
	case 'F':
	  if(argc<3)
	    {
	      printf("No file start given!\n");
	      usage();
	    }
	  strcpy(filestart, argv[2]);
	  argc--;
	  argv++;
	  break;
	case 'd':
	  if(argc<3)
	    {
	      printf("No smoothing time difference given!\n");
	      usage();
	    }
	  
	  dosmooth=1;
	  diff=atof(argv[2]);
	  argc--;
	  argv++;
	  break;
	case 'D':
	  if(argc<5)
	    {
	      printf("No smoothing time difference, start time, end time given!\n");
	      usage();
	    }
	  
	  dosmooth=1;
	  diff=atof(argv[2]);
	  t_start=atof(argv[3]);
	  t_end=atof(argv[4]);
	  argc-=3;
	  argv+=3;
	  break;
	case 'G':
	  {
	    if(argc<5)
	      {
		printf("No smoothing time difference, start time, end time given!\n");
		usage();
	      }
	    
	    dosmooth=1;
	    diff=atof(argv[2]);
	    HydDateTime dt1(argv[3]);
	    HydDateTime dt2(argv[4]);
	    dt_start=dt1;
	    dt_end=dt2;
	    do_dt_start_end=true;
	    
	    argc-=3;
	    argv+=3;
	    break;
	  }
	case 'R':
	  if(argc<3)
	    {
	      printf("No site to remove given!\n");
	      usage();
	    }
	  
	  use_site[atoi(argv[2])]=false;
	  
	  argc--;
	  argv++;
	  break;
	case 'O':
	  if(argc<3)
	    {
	      printf("No site to keep given!\n");
	      usage();
	    }
	  
	  for(i=0;i<1000;i++)
	    use_site[i]=false;
	  use_site[atoi(argv[2])]=true;
	  
	  argc--;
	  argv++;
	  break;
	case 't':
	  { 
	    switch(argv[1][2])
	      {
	      case 'l':
		talkative_likelihood=1;
		break;
	      case 'b': 
		talkative_burnin=1;
		break;
	      default:
		printf("Unknown talkative option!\n");
		usage();
		break;
	      }
	    break;
	  }
	case 'P':
	  { 
	    switch(argv[1][2])
	      {
	      case 'r':
		plot_repar=true;
		break;
	      case 's': 
		plot_sd=true;
		break;
	      default:
		printf("Unknown plotting option!\n");
		usage();
		break;
	      }
	    break;
	  }
	case 'C':
	  if(num_series_corr==0)
	    {
	      corr_from_series=new int[10000];
	      corr_to_series=new int[10000];
	      corr_from_layer=new int[10000];
	      corr_to_layer=new int[10000];
	      corr_from_index=new int[10000];
	      corr_to_index=new int[10000];
	      series_corr=new double[10000];
	    }
	  
	  corr_from_series[num_series_corr]=atoi(argv[2]);
	  corr_from_layer[num_series_corr]=atoi(argv[3]);
	  corr_to_series[num_series_corr]=atoi(argv[4]);
	  corr_to_layer[num_series_corr]=atoi(argv[5]);
	  
	  if(corr_from_series[num_series_corr]>(int)num_series)
	    {
	      cerr << "Series index too high:" << 
		corr_from_series[num_series_corr] << endl;
	      exit(0);
	    }

	  if(corr_to_series[num_series_corr]>(int)num_series)
	    {
	      cerr << "Series index too high:" << 
		corr_to_series[num_series_corr] << endl;
	      exit(0);
	    }

	  if(corr_from_series[num_series_corr]<=0)
	    {
	      cerr << "Series index too low:" << 
		corr_from_series[num_series_corr] << endl;
	      exit(0);
	    }

	  if(corr_to_series[num_series_corr]<=0)
	    {
	      cerr << "Series index too low:" << 
		corr_to_series[num_series_corr] << endl;
	      exit(0);
	    }

	  if(corr_from_layer[num_series_corr]>
	     (int)ser[corr_from_series[num_series_corr]-1].numlayers)
	    {
	      cerr << "Layer index too high:" << 
		corr_from_layer[num_series_corr] << endl;
	      exit(0);
	    }

	  if(corr_to_layer[num_series_corr]>
	     (int)ser[corr_to_series[num_series_corr]-1].numlayers)
	    {
	      cerr << "Layer index too high:" << 
		corr_to_layer[num_series_corr] << endl;
	      exit(0);
	    }

	  if(corr_from_layer[num_series_corr]<=0)
	    {
	      cerr << "Layer index too low:" << 
		corr_from_layer[num_series_corr] << endl;
	      exit(0);
	    }

	  if(corr_to_layer[num_series_corr]<=0)
	    {
	      cerr << "Layer index too low:" << 
		corr_to_layer[num_series_corr] << endl;
	      exit(0);
	    }
	  
	  corr_from_series[num_series_corr]--;
	  corr_from_layer[num_series_corr]--;
	  corr_to_series[num_series_corr]--;
	  corr_to_layer[num_series_corr]--;
	  
	  argc-=4;
	  argv+=4;
	  
	  num_series_corr++;
	  break;
	  /*
	    case 'r':
	    {
	    switch(argv[1][2])
	    {
	    case 'C':
	    if(num_series_corr==0)
	    {
	    corr_from_series=new int[10000];
	    corr_to_series=new int[10000];
	    corr_from_layer=new int[10000];
	    corr_to_layer=new int[10000];
	    corr_from_index=new int[10000];
	    corr_to_index=new int[10000];
	    corr_from_to_regional=new int[10000];
	    series_corr=new double[10000];
	    }
		
	    corr_from_series[num_series_corr]=atoi(argv[2]);
	    corr_from_layer[num_series_corr]=atoi(argv[3]);
	    corr_to_series[num_series_corr]=atoi(argv[4]);
	    corr_to_layer[num_series_corr]=atoi(argv[5]);
	    corr_from_to_regional[num_series_corr]=1;
		
	    if(corr_from_series[num_series_corr]>num_series)
	    {
	    cerr << "Series index too high:" << 
	    corr_from_series[num_series_corr] << endl;
	    exit(0);
	    }
		
	    if(corr_to_series[num_series_corr]>num_series)
	    {
	    cerr << "Series index too high:" << 
	    corr_to_series[num_series_corr] << endl;
	    exit(0);
	    }
		
	    if(corr_from_series[num_series_corr]<=0)
	    {
	    cerr << "Series index too low:" << 
	    corr_from_series[num_series_corr] << endl;
	    exit(0);
	    }
		
	    if(corr_to_series[num_series_corr]<=0)
	    {
	    cerr << "Series index too low:" << 
	    corr_to_series[num_series_corr] << endl;
	    exit(0);
	    }
		
	    if(corr_from_layer[num_series_corr]>
	    ser[corr_from_series[num_series_corr]-1].numlayers)
	    {
	    cerr << "Layer index too high:" << 
	    corr_from_layer[num_series_corr] << endl;
	    exit(0);
	    }
		
	    if(corr_to_layer[num_series_corr]>
	    ser[corr_to_series[num_series_corr]-1].numlayers)
	    {
	    cerr << "Layer index too high:" << 
	    corr_to_layer[num_series_corr] << endl;
	    exit(0);
	    }
		
	    if(corr_from_layer[num_series_corr]<=0)
	    {
	    cerr << "Layer index too low:" << 
	    corr_from_layer[num_series_corr] << endl;
	    exit(0);
	    }

	    if(corr_to_layer[num_series_corr]<=0)
	    {
	    cerr << "Layer index too low:" << 
	    corr_to_layer[num_series_corr] << endl;
	    exit(0);
	    }
	  
	    corr_from_series[num_series_corr]--;
	    corr_from_layer[num_series_corr]--;
	    corr_to_series[num_series_corr]--;
	    corr_to_layer[num_series_corr]--;
		
	    argc-=4;
	    argv+=4;
		
	    num_series_corr++;
	    break;
	    default:
	    printf("Unknown regional option!\n");
	    usage();
	    break;
	    }
	    break;
	    }
	  */
	case 'g':
	  if(num_series_feed==0)
	    {
	      feed_from_series=new int[10000];
	      feed_to_series=new int[10000];
	      feed_from_layer=new int[10000];
	      feed_to_layer=new int[10000];
	      beta_feed=new double[10000];
	      feed_symmetric=new int[10000];
	    }
	  
	  feed_from_series[num_series_feed]=atoi(argv[2]);
	  feed_from_layer[num_series_feed]=atoi(argv[3]);
	  feed_to_series[num_series_feed]=atoi(argv[4]);
	  feed_to_layer[num_series_feed]=atoi(argv[5]);
	  feed_symmetric[num_series_feed]=1;

	  if(feed_from_series[num_series_feed]>(int)num_series)
	    {
	      cerr << "Series index too high:" << 
		feed_from_series[num_series_feed] << endl;
	      exit(0);
	    }

	  if(feed_to_series[num_series_feed]>(int)num_series)
	    {
	      cerr << "Series index too high:" << 
		feed_to_series[num_series_feed] << endl;
	      exit(0);
	    }

	  if(feed_from_series[num_series_feed]<=0)
	    {
	      cerr << "Series index too low:" << 
		feed_from_series[num_series_feed] << endl;
	      exit(0);
	    }

	  if(feed_to_series[num_series_feed]<=0)
	    {
	      cerr << "Series index too low:" << 
		feed_to_series[num_series_feed] << endl;
	      exit(0);
	    }

	  if(feed_from_layer[num_series_feed]>
	     (int)ser[feed_from_series[num_series_feed]-1].numlayers)
	    {
	      cerr << "Layer index too high:" << 
		feed_from_layer[num_series_feed] << endl;
	      exit(0);
	    }

	  if(feed_to_layer[num_series_feed]>
	     (int)ser[feed_to_series[num_series_feed]-1].numlayers)
	    {
	      cerr << "Layer index too high:" << 
		feed_to_layer[num_series_feed] << endl;
	      exit(0);
	    }

	  if(feed_from_layer[num_series_feed]<=0)
	    {
	      cerr << "Layer index too low:" << 
		feed_from_layer[num_series_feed] << endl;
	      exit(0);
	    }

	  if(feed_to_layer[num_series_feed]<=0)
	    {
	      cerr << "Layer index too low:" << 
		feed_to_layer[num_series_feed] << endl;
	      exit(0);
	    }
	  
	  feed_from_series[num_series_feed]--;
	  feed_from_layer[num_series_feed]--;
	  feed_to_series[num_series_feed]--;
	  feed_to_layer[num_series_feed]--;
	  
	  argc-=4;
	  argv+=4;
	  
	  num_series_feed++;
	  break;
	case 'f':
	  if(num_series_feed==0)
	    {
	      feed_from_series=new int[10000];
	      feed_to_series=new int[10000];
	      feed_from_layer=new int[10000];
	      feed_to_layer=new int[10000];
	      beta_feed=new double[10000];
	      feed_symmetric=new int[10000];
	    }
	  
	  feed_from_series[num_series_feed]=atoi(argv[2]);
	  feed_from_layer[num_series_feed]=atoi(argv[3]);
	  feed_to_series[num_series_feed]=atoi(argv[4]);
	  feed_to_layer[num_series_feed]=atoi(argv[5]);
	  feed_symmetric[num_series_feed]=0;
	  
	  if(feed_from_series[num_series_feed]>(int)num_series)
	    {
	      cerr << "Series index too high:" << 
		feed_from_series[num_series_feed] << endl;
	      exit(0);
	    }

	  if(feed_to_series[num_series_feed]>(int)num_series)
	    {
	      cerr << "Series index too high:" << 
		feed_to_series[num_series_feed] << endl;
	      exit(0);
	    }

	  if(feed_from_series[num_series_feed]<=0)
	    {
	      cerr << "Series index too low:" << 
		feed_from_series[num_series_feed] << endl;
	      exit(0);
	    }

	  if(feed_to_series[num_series_feed]<=0)
	    {
	      cerr << "Series index too low:" << 
		feed_to_series[num_series_feed] << endl;
	      exit(0);
	    }

	  if(feed_from_layer[num_series_feed]>
	     (int)ser[feed_from_series[num_series_feed]-1].numlayers)
	    {
	      cerr << "Layer index too high:" << 
		feed_from_layer[num_series_feed] << endl;
	      exit(0);
	    }

	  if(feed_to_layer[num_series_feed]>
	     (int)ser[feed_to_series[num_series_feed]-1].numlayers)
	    {
	      cerr << "Layer index too high:" << 
		feed_to_layer[num_series_feed] << endl;
	      exit(0);
	    }

	  if(feed_from_layer[num_series_feed]<=0)
	    {
	      cerr << "Layer index too low:" << 
		feed_from_layer[num_series_feed] << endl;
	      exit(0);
	    }

	  if(feed_to_layer[num_series_feed]<=0)
	    {
	      cerr << "Layer index too low:" << 
		feed_to_layer[num_series_feed] << endl;
	      exit(0);
	    }
	  
	  feed_from_series[num_series_feed]--;
	  feed_from_layer[num_series_feed]--;
	  feed_to_series[num_series_feed]--;
	  feed_to_layer[num_series_feed]--;
	  
	  argc-=4;
	  argv+=4;
	  
	  num_series_feed++;
	  break;
	case 'c':
	  obs_corr=get_correlations(argv[2], &obs_corr_len, 
				    &obs_corr_site, &obs_corr_t, 
				    (ser[0].meas[0].dt!=NoHydDateTime ? 
				     true : false),
				    use_site, start, end);
	  argc--;
	  argv++;
	  break;
	case 'b':
	  do_importance=false;
	  break;
	case 'B':
	  if(argc<3)
	    {
	      printf("No start time given!\n");
	      usage();
	    }
	  
	  start=atof(argv[2]);
	  
	  argc--;
	  argv++;
	  break;
	case 'A':
	  if(argc<3)
	    {
	      printf("No end time given!\n");
	      usage();
	    }
	  
	  end=atof(argv[2]);
	  
	  argc--;
	  argv++;
	  break;
	case 's':
	  {
	    switch(argv[1][2])
	      {
	      case '\0':
		silent=1;
		break;
	      default:
		printf("Unknown seeparate option!\n");
		usage();
		break;
	      }
	    break;
	  }
	case 'e':
	  if(argc<3)
	    {
	      printf("No external time serie file given!\n");
	      usage();
	    }
	  
	  extdata=get_2d_data_file(argv[2], &ext_len);
	  if(ext_len<=0 || !extdata)
	    {
	      printf("Failed to read external time series from file '%s'!\n",
		     argv[2]);
	      exit(0);
	    }
	  useext=1;
	  
	  argc--;
	  argv++;
	  break;
	case 'E':
	  sel_pres=1;
	  break;
	case 'M':
	  if(argv[1][2]=='L')
	    do_ml=true;
	  else
	    {
	      printf("Unknown run option!\n");
	      usage();
	    }
	  
	  if(argc<3)
	    {
	      printf("No number of optimizations given!\n");
	      usage();
	    }
	  
	  num_optim=atoi(argv[2]);
	  
	  argc--;
	  argv++;
	  
	  break;
	default:
	  printf("Unknown run option!\n");
	  usage();
	  break;
	}
      
      argc--;
      argv++;
    }
  
  // if the wrong number of input, show usage:
  if(argc!=5 && !do_simulations)
    usage();
  
  num_states=0;
  num_tot_layers=0;
  for(s=0;s<num_series;s++)
    {
      for(l=0;l<ser[s].numlayers;l++)
	for(j=0;j<(int)numsites;j++)
	  {
	    state_series[num_states+l*numsites+j]=s;
	    state_layer[num_states+l*numsites+j]=l;
	  }
      
      series_state_start[s]=num_states;
      num_states += ser[s].numlayers*numsites;
      num_tot_layers+=ser[s].numlayers;
    }
  if(!silent)
    cout << "num_states=" << num_states << endl;

  if(num_series_corr>0)
    {
      for(i=0;i<num_series_corr;i++)
	{
	  int i1=series_state_start[corr_from_series[i]]+numsites*corr_from_layer[i];
	  int i2=series_state_start[corr_to_series[i]]+numsites*corr_to_layer[i]; 
	  
	  corr_from_index[i]=i1/numsites;
	  corr_to_index[i]=i2/numsites;
	}
    }

  if(num_series_corr>1)
    {
      // Estimate the probability of a correlation drawn from the prior
      // resulting in a positively definite covariance matrix.
      // Used for correcting this prior to this restriction.
      long int t3=clock();
      p_pos_series_sigma2=p_positive_seriescorr(10000);
      cout << "p_pos_series_sigma2=" << p_pos_series_sigma2 << endl;
      long int t4=clock();
      cout << "Computing time=" << double(t4-t3)/
	double(CLOCKS_PER_SEC) << endl;
      //exit(0);
    }
  else
    p_pos_series_sigma2=1.0;


  meas_tot_len=0;
  for(s=0;s<num_series;s++)
    {
      meas_tot_len += ser[s].meas_len;
      for(i=0;i<ser[s].meas_len;i++)
	ser[s].meas[i].index=series_state_start[s]+ser[s].meas[i].site;
    }
  series_measurements *meas_buffer=new series_measurements[meas_tot_len];	  
  j=0;
  for(s=0;s<num_series;s++)
    for(i=0;i<ser[s].meas_len;i++)
      meas_buffer[j++]=ser[s].meas[i];
  qsort(meas_buffer,(size_t) meas_tot_len,sizeof(series_measurements), 
	compare_meas);
  
  meas_tot=new measurement_cluster[meas_tot_len];
  j=-1;
  for(i=0;i<meas_tot_len;i++)
    {
      if(j<0 || meas_tot[j].tm!=meas_buffer[i].tm)
	j++;
      if(meas_buffer[i].meanval!=MISSING_VALUE)
	meas_tot[j].add_measurement(meas_buffer[i]);
      else if(meas_tot[j].num_measurements==0)
	{
	  meas_tot[j].tm=meas_buffer[i].tm;
	  meas_tot[j].dt=meas_buffer[i].dt;
	}
      //if(meas_tot[j].num_measurements>1)
      //meas_tot[j].print();
      
      //ser[s].meas[i].print();
      //meas_tot[j].print();
    }
  delete [] meas_buffer;
  meas_tot_len=j+1;
  printf("meas_tot_len=%d\n", meas_tot_len);
  //qsort(meas_tot,(size_t) meas_tot_len,sizeof(measurement_cluster), 
  //compare_meas_cluster);
  if(do_dt_start_end)
    {
      int secs1=dt_start-meas_tot[0].dt;
      int secs2=dt_end-meas_tot[0].dt;
      t_start=(double)secs1;
      t_end=(double)secs2;
    }


  // find the minimal time difference:
  double min_time_diff=meas_tot[meas_tot_len-1].tm-meas_tot[0].tm;
  for(i=1;i<meas_tot_len;i++)
    min_time_diff=MINIM(min_time_diff, (meas_tot[i].tm-meas_tot[i-1].tm));
  if(min_time_diff<=0.0)
    {
      cerr << "min_time_diff=" << min_time_diff << "!" << endl;
      exit(0);
    }

  // If there is correlation information, stitch that in:
  for(i=0;i<meas_tot_len;i++)
    {
      int *handled=new int[meas_tot[i].num_measurements];
      for(j=0;j<(int)meas_tot[i].num_measurements;j++)
	handled[j]=0;
      
      double **total_corr_matrix=Make_matrix(meas_tot[i].num_measurements,
					     meas_tot[i].num_measurements);
      
      for(j=0;j<(int)meas_tot[i].num_measurements;j++)
	if(!handled[j])
	  {
	    // site+time associated with observational correlation matrix?
	    
	    for(k=0;(int)k<obs_corr_len;k++)
	      if(obs_corr_t[k]>=(meas_tot[i].tm-min_time_diff/2) &&
		 obs_corr_t[k]<=(meas_tot[i].tm+min_time_diff/2))
		break;
	    if((int)k<obs_corr_len) // observational correlation matrix found
	      {
		// Make observational correlation matrix for this site+time:
		double **corr_matrix=Make_matrix(num_series,num_series);
		for(l=0;l<num_series;l++)
		  corr_matrix[l][l]=1.0;
		m=0;
		for(l=0;l<(num_series-1);l++)
		  for(s=l+1;s<num_series;s++)
		    {
		      corr_matrix[l][s]=corr_matrix[s][l]=obs_corr[k][m];
		      m++;
		    }
		
		for(k=0;k<meas_tot[i].num_measurements;k++)
		  for(l=0;l<meas_tot[i].num_measurements;l++)
		    if(meas_tot[i].site[k]==meas_tot[i].site[j] &&
		       meas_tot[i].site[l]==meas_tot[i].site[j])
		      {
			total_corr_matrix[k][l]=
			  corr_matrix[meas_tot[i].serie_num[k]][meas_tot[i].serie_num[l]];
			handled[k]=handled[l]=1;
		      }

		doubledelete(corr_matrix, num_series);
	      }
	  }
      
      for(j=0;j<(int)meas_tot[i].num_measurements;j++)
	total_corr_matrix[j][j]=1.0;
      
      meas_tot[i].corr_matrix=Make_matrix(meas_tot[i].num_measurements,
					  meas_tot[i].num_measurements);
      for(j=0;j<(int)meas_tot[i].num_measurements;j++)
	for(k=0;k<meas_tot[i].num_measurements;k++)
	  meas_tot[i].corr_matrix[j][k]=total_corr_matrix[j][k];

      doubledelete(total_corr_matrix,meas_tot[i].num_measurements);
    }


  //if(!silent)
  //for(i=0;i<meas_tot_len;i++)
  //  meas_tot[i].print();

  HydDateTime ref(1970);
  /*
    for(j=0;j<meas_tot_len;j++)
    {
    HydDateTime dt=ref+((long int)meas_tot[j].tm);
    printf("%f %s %s\n", meas_tot[j].tm, meas_tot[j].dt.syCh(1),
    dt.syCh(1));
    }
  */

  // Make measurement clusters for smoothed time series:
  if(dosmooth)
    {
      if(diff<=0.0)
	{
	  meas_smooth_len=meas_tot_len;
	  meas_smooth=new measurement_cluster[meas_tot_len];
	  for(i=0;i<meas_tot_len;i++)
	    meas_smooth[i].copy(&(meas_tot[i]));
	}
      else
	{
	  meas_smooth_len=0;
	  double t1=(t_start==MISSING_VALUE ? meas_tot[0].tm-1.0 : t_start);
	  double t2=(t_end==MISSING_VALUE ? 
		     MAXIM(0.0,meas_tot[meas_tot_len-1].tm+1.0) :
		     t_end);
	  double t=t1;
      
	  if(meas_tot[0].dt!=NoHydDateTime)
	    {
	      if(t_start==MISSING_VALUE)
		t1=meas_tot[0].tm-
		  0.1*double(meas_tot[meas_tot_len-1].dt.
			     difference_minutes(meas_tot[0].dt))*60.0;

	      if(t_end==MISSING_VALUE)
		t2=meas_tot[meas_tot_len-1].tm+
		  0.1*double(meas_tot[meas_tot_len-1].dt.
			     difference_minutes(meas_tot[0].dt))*60.0;
	      t=t1;
	    }
	  
	  for(t=t1;t<=t2;t+=diff)
	    meas_smooth_len++;
	  cout << "Meas_smooth_len=" << meas_smooth_len << endl;
	  
	  meas_smooth=new measurement_cluster[meas_smooth_len+meas_tot_len+10];
	  unsigned int j=0,k=0;

	  i=0;
	  t=MINIM(t1,meas_tot[0].tm);
	  while(i<meas_smooth_len+meas_tot_len+1 && 
		(t<MAXIM(t2,meas_tot[meas_tot_len-1].tm)))
	    {
	      while(j<meas_tot_len && meas_tot[j].tm<t)
		{
		  //cout << k << " av " << meas_smooth_len+meas_tot_len+10 << " "
		  //   << j << " av " << meas_tot_len << endl;
		  
		  meas_smooth[k].copy(&(meas_tot[j]));
		  //cout << k << " " << meas_smooth[k].dt << " " << 
		  //meas_tot[j].dt << " " << meas_tot[0].dt << endl;
		  k++;
		  j++;
		}
	      
	      if((j>=meas_tot_len || meas_tot[j].tm>t) &&
		 (t>=t1 && t<=t2))
		{
		  measurement_cluster buff(t);
		  meas_smooth[k].copy(&buff);
		  if(meas_tot[0].dt!=NoHydDateTime)
		    {
		      meas_smooth[k].dt=ref+((long int)t);
		      //cout << ref << " " << t << " " << ((long int) t) << endl;
		      //cout << k << " " << meas_smooth[k].dt << " " << 
		      //meas_tot[0].dt << endl;
		      if(!meas_smooth[k].dt.legal())
			{
			  cout << ref << " " << t << " " << ((long int)t) << 
			    " " << meas_tot[0].dt << endl;
			  exit(0);
			}
		    }
		  k++;
		}

	      i++;
	      t=t1+diff*double(i);
	    }
	  
	  meas_smooth_len=k;
	  cout << meas_smooth_len << endl;
	}
      
      /*
	if(!silent)
	for(i=0;i<meas_smooth_len;i++)
	meas_smooth[i].print();
      */

      //for(i=0;i<meas_smooth_len;i++)
      //printf("%s\n",meas_smooth[i].dt.syCh(1));
    }
  
  


  //for(i=0;i<meas_smooth_len;i++)
  //meas_smooth[i].print();

  if(do_simulations)
    {
      if(argc!=1)
	usage();
      
      double *sim_param_repar=new double[numpar];
      for(i=0;i<numpar;i++)
	sim_param_repar[i]=transform_parameter(sim_param[i],par_trans_type[i]);
      
      char **filenames=new char*[num_series];
      for(s=0;s<num_series;s++)
	filenames[s]=new char[1000];
      
      double ll;
      for(i=0;(int)i<numsim;i++)
	{
	  if(num_series>0)
	    for(s=0;s<num_series;s++)
	      snprintf(filenames[s], 800,"%s_serie%02d_sim%04d.txt", 
		      simulation_file_start,s,i+1);
	  else
	    snprintf(filenames[s], 800,"%s_sim%04d.txt",
		     simulation_file_start,i+1);
	  
	  ll=loglik(sim_param_repar, 0, 0, 0, NULL, 0, filenames);
	}
      cout << "loglik=" << ll << endl;

      doubledelete(filenames,num_series);
      exit(0);
    }


  if(num_series<1)
    {
      cerr << "No series given!" << endl;
      exit(0);
    }
  
  if(do_realizations && !dosmooth)
    {
      printf("Realizations sent to file specified but no smoothing!\n");
      exit(0);
    }
  
  // read the number of parameters and also fill out the global
  // variables par_name and par_type:
  loglik(NULL);
  
  // fetch sampling parameters:
  int numsamples=atoi(argv[1]);
  int burnin=atoi(argv[2]);
  int indep=atoi(argv[3]);
  unsigned int numtemp=(unsigned int)atoi(argv[4]);
    
  if(preshow)
    for(s=0;s<num_series;s++)
      for(j=0;j<(int)numsites;j++)
	{
	  char cmd[1000];
	  
	  if(ser[s].meas[0].dt!=NoHydDateTime)
	    snprintf(cmd,999,"timeseriegraph -x t -y %s",ser[s].name);
	  else
	    snprintf(cmd,999,"vvgraph -x t -y %s",ser[s].name);
	  
	  FILE *p;
	  p=popen(cmd,"w");
	  fprintf(p, "# Column 1: measurements\n");
	  fprintf(p, "# Column 1 - type: dot\n");
	  if(ser[s].comp_meas)
	    {
	      fprintf(p, "# Column 2: comparison measurements\n");
	      fprintf(p, "# Column 2 - type: dot\n");
	    }
	  fprintf(p, "########################\n");
	  s=state_series[i];
	  for(k=0;k<ser[s].meas_len;k++)
	    if(ser[s].meas[k].site==(int)j)
	      {
		if(ser[s].meas[0].dt!=NoHydDateTime)
		  {
		    if(!ser[s].pr->is_log)
		      fprintf(p,"%s %f", ser[s].meas[k].dt.syCh(5), 
			      ser[s].meas[k].meanval);
		    else
		      fprintf(p,"%s %f", ser[s].meas[k].dt.syCh(5), 
			      ser[s].meas[k].meanval!=MISSING_VALUE ?
			      exp(ser[s].meas[k].meanval) : MISSING_VALUE);
		    
		  }
		else
		  {
		    if(!ser[s].pr->is_log)
		      fprintf(p,"%f %f", ser[s].meas[k].tm, ser[s].meas[k].meanval);
		    else
		      fprintf(p,"%f %f", ser[s].meas[k].tm, 
			      ser[s].meas[k].meanval!=MISSING_VALUE ?
			      exp(ser[s].meas[k].meanval) : MISSING_VALUE);
		    
		  }
	      
		if(ser[s].comp_meas)
		  fprintf(p," -10000000");
		fprintf(p,"\n");
	      }
	  
	  if(ser[s].comp_meas)
	    for(k=0;k<ser[s].comp_len;k++)
	      if(ser[s].comp_meas[k].site==(int)j)
		{
		  if(ser[s].meas[0].dt!=NoHydDateTime)
		    {
		      if(!ser[s].pr->is_log)
			fprintf(p,"%s -10000000 %f\n", ser[s].comp_meas[k].dt.syCh(5), 
				ser[s].comp_meas[k].meanval);
		      else
			fprintf(p,"%s -10000000 %f\n", ser[s].meas[k].dt.syCh(5), 
				ser[s].comp_meas[k].meanval!=MISSING_VALUE ?
				exp(ser[s].comp_meas[k].meanval) : MISSING_VALUE);
		    }
		  else
		    {
		      if(!ser[s].pr->is_log)
			fprintf(p,"%f -10000000 %f\n", 
				ser[s].comp_meas[k].tm, ser[s].meas[k].meanval);
		      else
			fprintf(p,"%f -10000000 %f\n", ser[s].meas[k].tm, 
				ser[s].comp_meas[k].meanval!=MISSING_VALUE ?
				exp(ser[s].comp_meas[k].meanval) : MISSING_VALUE);
		    }
		}
	  pclose(p);
	}
  
  
  /*
  // DEBUG for a specific model
  double par1[]={1.990356, log(0.155857), log(1.315038), log(0.105345),
  2.0, 2.015047};
  double par2[]={1.990356, log(0.155857), log(1.315038), log(0.105345),
  -2655392692221411860022296576.000000, 2.015047};
  params p1(par1,14), p2(par2,14);

  cout << logprob(p1, 1.0) << " " << loglik(par1,0,0,1) << endl;
  cout << logprob(p2, 1.0) << " " << loglik(par2,0,0,1) << endl;

  exit(0);
  */

  /*
  // DEBUG for a specific model
  double par1[]={log(7.426756), log(0.000527), log(0.305514),
  log(0.004881), log(0.677918), log(1.140461),
  log(0.029075), log(1.042266), log(1.231357),
  log(0.020987), log(3.753786), log(0.136665),
  transform_parameter(0.488898, T_LOGIST_GLOBAL),
  0.0};
  
  double par2[]={log(7.426756), log(0.000527), log(0.305514),
  log(0.004881), log(0.677918), log(1.140461),
  log(0.029075), log(1.042266), log(1.231357),
  log(0.020987), log(3.753786), log(0.136665),
  transform_parameter(0.488898, T_LOGIST_GLOBAL),
  1.0};
  params p1(par1,14), p2(par2,14);
  
  cout << logprob(p1, 1.0) << " " << loglik(par1) << endl;
  cout << logprob(p2, 1.0) << " " << loglik(par2) << endl;
  */
  
  
  double ***x; // smoothing results from the latent processes
  
  // fetch mcmc samples (also calculates the marginal
  // data density):
  long int t1=clock();
  double model_loglik=MISSING_VALUE;
  params *pars=layer_mcmc(numsamples, burnin, indep, numtemp, 
			  do_importance, dosmooth, do_realizations, 
			  realization_file_start,&x,startpar, T_ground, 
			  &model_loglik);
  long int t2=clock();
  
  cout << "Computing time=" << double(t2-t1)/
    double(CLOCKS_PER_SEC) << "s" << endl;
  
  //cout << "T0=" << T0 << " T1=" << T1 << " T2=" << T2 << " T3=" << T3 << 
  //" T4=" << T4 << " T5=" << T5 << " T6=" << T6 << endl;
  
  // *********************************
  // Show the results:
  
  // Inverse transform the parameter samples back to
  // original parameterization:
  double **parsample=Make_matrix(numpar,numsamples);
  double **parsample_repar=Make_matrix(numpar,numsamples);
  for(i=0;i<numpar;i++)
    for(j=0;(int)j<numsamples;j++)
      {
	parsample[i][j]=invtransform_parameter(pars[j].param[i], par_trans_type[i]);
	parsample_repar[i][j]=pars[j].param[i];
      }
      
  char **par_name_orig=new char*[numpar];
  for(i=0;i<numpar;i++)
    {
      par_name_orig[i]=new char[300];
      strcpy(par_name_orig[i],par_name[i]);
    }
  
  // transform the expectancies from logarithmic size to
  // original size:
  for(i=0;i<numpar;i++)
    {
      s=par_series[i];
      if(ser[s].pr->is_log)
	{
	  if(par_type[i]==MU)
	    {
	      snprintf(par_name[i],99,"exp(%s)",par_name_orig[i]);
	      for(j=0;(int)j<numsamples;j++)
		parsample[i][j]=exp(parsample[i][j]);
	    }
	  if(par_type[i]==OBS_SD)
	    {
	      snprintf(par_name[i],99,"%s_origscale",par_name_orig[i]);
	      for(j=0;(int)j<numsamples;j++)
		parsample[i][j]=(ser[s].pr->is_log==2 ? 
				 exp(ser[s].mean_val) : ser[s].mean_val)*
		  sqrt(exp(parsample[i][j]*parsample[i][j])-1.0)*
		  exp(parsample[i][j]*parsample[i][j]/2.0);
	    }
	}
    }


  
  // show mean, median and 95% credibility interval for
  // each parameter:
  for(j=0;j<(int)numpar;j++)
  {
    double mean=find_statistics(parsample[j], numsamples, MEAN);
    double med=find_statistics(parsample[j], numsamples, MEDIAN);
    double lower=find_statistics(parsample[j], numsamples, PERCENTILE_2_5);
    double upper=find_statistics(parsample[j], numsamples, PERCENTILE_97_5);
    
    if(report_halflife && !strncmp(par_name[j],"dt_",3))
      {
	mean*=log(2.0);
	med*=log(2.0);
	lower*=log(2.0);
	upper*=log(2.0);
      }
    
    printf("%s: mean=%f median=%f    95%% cred=(%f-%f)\n", 
	   par_name[j], mean, med, lower,upper);
  }
  
 
  if(do_ml)
    {
      // traverse the number of wanted hill-climbs:
      double *best_pars_repar=new double[numpar];
      double *best_pars=new double[numpar];
      double best_loglik=MISSING_VALUE;
      
      for(i=0;(int)i<num_optim;i++)
	{
	  double *curr_par=new double[numpar];
	  
	  for(j=0;j<(int)numpar;j++)
	    {
	      if(i==0)
		curr_par[j]=find_statistics(parsample_repar[j],numsamples,MEDIAN);
	      else if(i==1)
		curr_par[j]=find_statistics(parsample_repar[j],numsamples,MEAN);
	      else if(num_optim==3)
		curr_par[j]=parsample_repar[j][numsamples/2];
	      else
		curr_par[j]=parsample_repar[j][(i-2)*(numsamples-1)/(num_optim-3)];
	    }
	  
	  // do the optimization:
	  int maxiter=1000;
	  //double *pars2=gsl_optimization_cover(minusloglik, 
	  double *pars2=quasi_newton(minusloglik,
				     numpar /* number of parameters */, 
				     curr_par /* starting values */, 
				     0.001 /* precision */, 
				     maxiter /* Maximum number of
						iterations. The number
						of iterations
						neccessary will be
						stored here after the
						optimization and can 
						be checked */);
	  
	  // Fetch the log-likelihood for the optimized parameters:
	  double curr_loglik=loglik(pars2);
	  
	  // Check if this is the best result so far:
	  if(i==0 || curr_loglik>best_loglik)
	    {
	      // If so, store the parameter array and the log-likelihood:
	      for(j=0;j<(int)numpar;j++)
		best_pars_repar[j]=pars2[j];
	      best_loglik=curr_loglik;
	    }
	  
	  delete [] curr_par;
	  delete [] pars2;
	}
      
      // transform the expectancies from logarithmic size ot
      // original size:
      for(j=0;j<(int)numpar;j++)
	{
	  best_pars[j]=invtransform_parameter(best_pars_repar[j], 
					      par_trans_type[j]);
	  
	  s=par_series[j];
	  if(ser[s].pr->is_log)
	    {
	      if(par_type[j]==MU)
		best_pars[j]=exp(best_pars[j]);
	      if(par_type[j]==OBS_SD)
		best_pars[j]=ser[s].mean_val*
		  sqrt(exp(best_pars[j]*best_pars[j])-1.0)*
		  exp(best_pars[j]*best_pars[j]/2.0);
	    }
	}
      
      printf("\nML: %g\n", best_loglik);
      for(j=0;j<(int)numpar;j++)
	printf("%s - ML: %f \n", 
	       par_name[j], best_pars[j]);
      printf("\n");
      
      loglik(best_pars_repar, 0, 0,1,filestart);
      
      delete [] best_pars;
      delete [] best_pars_repar;
      
      k=numpar;
      double n=0.0;
      for(s=0;s<num_series;s++)
	n+=double(ser[s].meas_len);
      printf("-AIC/2  = %f\n", best_loglik - double(k));
      printf("-AICc/2 = %f\n", best_loglik - double(k) -
	     double(k*(k+1))/double(n-k-1)); 
      printf("-BIC/2  = %f\n", best_loglik - 0.5*double(k)*log(double(n)));
      printf("%f\n",log(double(n)));
    }

  if(num_series_feed>0)
    {
      double *iscomplex=new double[numsamples];

      for(j=0;(int)j<numsamples;j++)
	iscomplex[j]=pars[j].iscomplex ? 1.0 : 0.0;

      printf("%s: mean=%f median=%f    95%% cred=(%f-%f)\n", 
	     "iscomplex", 
	     find_statistics(iscomplex, numsamples, MEAN),
	     find_statistics(iscomplex, numsamples, MEDIAN),
	     find_statistics(iscomplex, numsamples, PERCENTILE_2_5),
	     find_statistics(iscomplex, numsamples, PERCENTILE_97_5));
      show_parameter(iscomplex, numsamples, (char *) "iscomplex", silent, filestart);

      delete [] iscomplex;

      for(i=0;i<MINIM(ser[0].numlayers,10);i++)
	{
	  double *cycle=new double[numsamples];
	  for(j=0;(int)j<numsamples;j++)
	    cycle[j]=pars[j].cycles[i];

	  printf("%s%02d: mean=%f median=%f    95%% cred=(%f-%f)\n", 
		 "cycle",i+1, 
		 find_statistics(cycle, numsamples, MEAN,0.01),
		 find_statistics(cycle, numsamples, MEDIAN,0.01),
		 find_statistics(cycle, numsamples, PERCENTILE_2_5,0.01),
		 find_statistics(cycle, numsamples, PERCENTILE_97_5,0.01));
	  char name[100];
	  snprintf(name,99,"cycle%02d",i+1);
	  show_parameter(cycle, numsamples, name, silent, filestart);
	  
	  delete [] cycle;
	}
    }
  
  // if not in silent modus, show inferred series:
  if(!silent && dosmooth && !do_realizations)
    {
      char filename[1000], cmd[2000];
      
      if(ser[0].numlayers>1 && ser[0].meas[0].sd!=MISSING_VALUE && sel_pres)
	// check measurement variance vs selection pressure (difference
	// between layer one and two) for the base series
	{
	  FILE *p=popen("vvgraph -x \"selection pressure\" -y sd","w");
	  fprintf(p,"# Column 1: measurement standard deviation "
		  "vs selection pressure\n");
	  fprintf(p,"# Column 1 - type : dot\n");
	  fprintf(p,"###########################\n");
	  for(k=0;k<meas_smooth_len;k++)
	    if(meas_smooth[k].num_measurements>0)
	      {
		for(m=0;m<meas_smooth[k].num_measurements;m++)
		  {
		    if(meas_smooth[k].serie_num[m]==0 && 
		       meas_smooth[k].meanval[m]!=MISSING_VALUE)
     {
			// site:
			j=meas_smooth[k].site[m];
			
			double upper_layer=
			  find_statistics(x[j][k], 
					  num_smooth*numsamples, MEAN);
			double lower_layer=
			  find_statistics(x[numsites+j][k], 
					  num_smooth*numsamples, MEAN);
			
			double selection_pressure=ABSVAL((upper_layer-lower_layer));
			
			fprintf(p, "%f %f\n", selection_pressure, meas_smooth[k].sd[0]); 
		      }
		  }
	      }
	  pclose(p);
	}
      
      // Layer 0, all sites
      if(numsites>1)
	for(s=0;s<num_series;s++)
	  {
	    i=series_state_start[s];
	    l=0;
	    
	    if(ser[s].meas[0].dt!=NoHydDateTime)
	      snprintf(cmd,1500,"timeseriegraph -x t -y %s",ser[s].name);
	    else
	      snprintf(cmd,1500,"vvgraph -x t -y %s",ser[s].name);
	    
	    FILE *p;
	    if(!*filestart)
	      p=popen(cmd,"w");
	    else
	      {
		snprintf(filename, 800, "%s_allsites_%s_layer0.txt", 
			filestart,ser[s].name);
		p=fopen(filename,"w");
	      }
	    for(j=0;j<(int)numsites;j++)
	      {
		fprintf(p, "# Column %d: mean %s, site %d, layer %d, x%d\n", 
			2*j+1,ser[s].name,j, l+1, j);
		fprintf(p, "# Column %d: Measurements of %s for site %d\n", 
			2*j+2,ser[s].name,j);
	      }
	    fprintf(p, "########################\n");
	    
	    for(k=0;k<meas_smooth_len;k++)
	      {
		/*
		  if(t_k_smooth[k]> -10.045 && t_k_smooth[k]< -10.035)
		  cout << find_statistics(x[0][k], 
		  num_smooth*numsamples, MEAN) << " " <<
		  find_statistics(x[1][k], 
		  num_smooth*numsamples, MEAN) << endl;
		*/
		
		if(ser[s].meas[0].dt!=NoHydDateTime)
		  fprintf(p,"%s ", dt_k_smooth[k].syCh(5));
		else
		  fprintf(p,"%f ", t_k_smooth[k]);
		
		for(j=0;j<(int)numsites;j++)
		  {
		    bool done=false;
		    for(m=0;m<=meas_smooth[k].num_measurements && !done;m++)
		      if(m==meas_smooth[k].num_measurements ||
			 (meas_smooth[k].meanval[m]!=MISSING_VALUE &&
			  meas_smooth[k].site[m]==(int)j && 
			  meas_smooth[k].serie_num[m]==(int)s))
			{
			  if(!ser[s].pr->is_log)
			    fprintf(p,"%f %f ",
				    find_statistics(x[i+j][k], 
						    num_smooth*numsamples, MEAN),
				    m<meas_smooth[k].num_measurements ? 
				    meas_smooth[k].meanval[m] : MISSING_VALUE);
			  else
			    fprintf(p,"%f %f ",
				    exp(find_statistics(x[i+j][k], 
							num_smooth*numsamples, MEAN)),
				    (m<meas_smooth[k].num_measurements) ? 
				    exp(meas_smooth[k].meanval[m]) : MISSING_VALUE);
			  done=true;
			}
		  }
		fprintf(p,"\n");
	      }
	    
	    if(!*filestart)
	      pclose(p);
	    else
	      fclose(p);
	  }
      
      
      
      // site+layer with uncertainty:
      for(i=0;i<num_states;i++)
	{
	  s=state_series[i];
	  l=state_layer[i];
	  
	  if(ser[s].meas[0].dt!=NoHydDateTime)
	    snprintf(cmd,999,"timeseriegraph -x t -y %s",ser[s].name);
	  else
	    snprintf(cmd,999,"vvgraph -x t -y %s",ser[s].name);
	  
	  FILE *p;
	  if(!*filestart)
	    p=popen(cmd,"w");
	  else
	    {
	      snprintf(filename,800, "%s_singlesite%02d.txt", filestart, i);
	      p=fopen(filename,"w");
	    }
	  fprintf(p, "# Column 1: mean %s, site %d, layer %d, x%d\n", 
		  ser[s].name, i%numsites, l+1, i);
	  fprintf(p, "# Column 2: lower 95%% cred\n");
	  fprintf(p, "# Column 3: lower 95%% cred\n");
	  fprintf(p, "# Column 4: measurements\n");
	  fprintf(p, "# Column 4 - type: dot\n");
	  if(ser[s].comp_meas)
	    {
	      fprintf(p, "# Column 5: comparison measurements\n");
	      fprintf(p, "# Column 5 - type: dot\n");
	    }
	  fprintf(p, "########################\n");
	  for(k=0;k<meas_smooth_len;k++)
	    {
	      bool done=false;
	      for(m=0;m<=meas_smooth[k].num_measurements && !done;m++)
		if(m==meas_smooth[k].num_measurements ||
		   (meas_smooth[k].meanval[m]!=MISSING_VALUE &&
		    meas_smooth[k].site[m]==(int)(i%numsites) &&
		    meas_smooth[k].serie_num[m]==(int)s))
		  {
		    if(ser[s].meas[0].dt!=NoHydDateTime)
		      {
			if(dt_k_smooth[k].legal())
			  fprintf(p,"%s ",dt_k_smooth[k].syCh(5));
			else
			  {
			    cerr << "Illegal time for time point " << 
			      k << " : " << dt_k_smooth[k] << "=" << 
			      t_k_smooth[k] << endl;
			    exit(0);
			  }
		      }
		    else
		      fprintf(p,"%f ",t_k_smooth[k]);
		    
		    if(!ser[s].pr->is_log)
		      fprintf(p,"%f %f %f %f", 
			      find_statistics(x[i][k], num_smooth*numsamples, MEAN),
			      find_statistics(x[i][k], num_smooth*numsamples, 
					      PERCENTILE_2_5),
			      find_statistics(x[i][k], num_smooth*numsamples, 
					      PERCENTILE_97_5),
			      (m<meas_smooth[k].num_measurements) ? 
			      meas_smooth[k].meanval[m] : MISSING_VALUE);
		    else
		      fprintf(p,"%f %f %f %f", 
			      exp(find_statistics(x[i][k], num_smooth*numsamples, MEAN)),
			      exp(find_statistics(x[i][k], num_smooth*numsamples, 
						  PERCENTILE_2_5)),
			      exp(find_statistics(x[i][k], num_smooth*numsamples, 
						  PERCENTILE_97_5)),
			      (m<meas_smooth[k].num_measurements) ? 
			      exp(meas_smooth[k].meanval[m]) : MISSING_VALUE);
		    
		    done=true;
		  }
	      if(ser[s].comp_meas)
		fprintf(p," -10000000");
	      fprintf(p,"\n");
	    }
	  
	  if(ser[s].comp_meas)
	    for(k=0;k<ser[s].comp_len;k++)
	      if(ser[s].comp_meas[k].site==(int)(i%numsites))
		{
		  if(ser[s].meas[0].dt!=NoHydDateTime)
		    {
		      if(!ser[s].pr->is_log)
			fprintf(p,"%s -10000000 -10000000 -10000000 -10000000 %f\n", 
				ser[s].comp_meas[k].dt.syCh(5), 
				ser[s].comp_meas[k].meanval);
		      else
			fprintf(p,"%s -10000000 -10000000 -10000000 -10000000 %f\n", 
				ser[s].meas[k].dt.syCh(5), 
				ser[s].comp_meas[k].meanval!=MISSING_VALUE ?
				exp(ser[s].comp_meas[k].meanval) : 
				MISSING_VALUE);
		    }
		  else
		    {
		      if(!ser[s].pr->is_log)
			fprintf(p,"%f -10000000 -10000000 -10000000 -10000000 %f\n", 
				ser[s].comp_meas[k].tm, ser[s].meas[k].meanval);
		      else
			fprintf(p,"%f -10000000 -10000000 -10000000 -10000000 %f\n", 
				ser[s].meas[k].tm, 
				ser[s].comp_meas[k].meanval!=MISSING_VALUE ?
				exp(ser[s].comp_meas[k].meanval) : 
				MISSING_VALUE);
		    }
		}
	  
	  if(!*filestart)
	    pclose(p);
	  else
	    fclose(p);
	}

      if(numsites>1)
	{ 
	  // layer by layer plot for all sites
	  for(s=0;s<num_series;s++)
	    for(l=0;l<ser[s].numlayers;l++)
	      {
		FILE *p;
		unsigned int j,k;
		
		if(ser[s].meas[0].dt!=NoHydDateTime)
		  snprintf(cmd,999,"timeseriegraph -x t -y %s",ser[s].name);
		else
		  snprintf(cmd,999,"vvgraph -x t -y %s",ser[s].name);
		
		if(!*filestart)
		  p=popen(cmd,"w");
		else
		  {
		    snprintf(filename,900, "%s_%s_layer%01d.txt", filestart, 
			    ser[s].name,l+1);
		    p=fopen(filename,"w");
		  }
		for(j=0;j<numsites;j++)
		  fprintf(p, "# Column %d: mean %s, x%d, site %d, layer %d\n", 
			  j+1, ser[s].name, 
			  series_state_start[s]+l*numsites+j,j,l+1);
		fprintf(p, "########################\n");
		for(k=0;k<meas_smooth_len;k++)
		  {
		    if(ser[s].meas[0].dt!=NoHydDateTime)
		      fprintf(p, "%s ",dt_k_smooth[k].syCh(5));
		    else
		      fprintf(p, "%f ",t_k_smooth[k]);
		    for(j=0;j<numsites;j++)
		      {
			if(!ser[s].pr->is_log)
			  fprintf(p, "%f ",
				  find_statistics(x[series_state_start[s]+
						    l*numsites+j][k], 
						  num_smooth*numsamples, MEAN));
			else
			  fprintf(p, "%f ",
				  exp(find_statistics(x[series_state_start[s]+
							l*numsites+j][k],
						      num_smooth*numsamples, MEAN)));
		      }
		    fprintf(p, "\n");
		  }
		
		if(!*filestart)
		  pclose(p);
		else
		  fclose(p);
	      }
	}
      


      // site by site showing all layers:
      for(s=0;s<num_series;s++)
	if(ser[s].numlayers>=2)
	  {
	    for(j=0;j<(int)numsites;j++)
	      {
		FILE *p;
		if(!*filestart)
		  {
		    if(ser[s].meas[0].dt!=NoHydDateTime)
		      p=popen("timeseriegraph", "w");
		    else
		      p=popen("vvgraph", "w");
		  }
		else
		  {
		    snprintf(filename, 800,"%s_site%01d.txt", filestart, j);
		    p=fopen(filename, "w");
		  }
		k=1;
		for(l=0;l<ser[s].numlayers;l++)
		  fprintf(p, "# Column %d: %s site %d, mean x%d\n", k++,
			  ser[s].name,j,
			  series_state_start[s]+l*numsites+j);
		fprintf(p,"# Column %d: measurements\n", k);
		fprintf(p,"# Column %d -type: dot\n", k);
		if(ser[s].comp_meas)
		  {
		    k++;
		    fprintf(p,"# Column %d: comparison measurements\n", k);
		    fprintf(p,"# Column %d -type: dot\n", k);
		  }
		  
		fprintf(p, "########################\n");
		for(k=0;k<meas_smooth_len;k++)
		  {
		    if(ser[s].meas[0].dt!=NoHydDateTime)
		      fprintf(p, "%s ",dt_k_smooth[k].syCh(5));
		    else
		      fprintf(p, "%f ",t_k_smooth[k]);

		    for(l=0;l<ser[s].numlayers;l++)
		      if(!ser[s].pr->is_log)
			fprintf(p, "%f ", 
				find_statistics(x[series_state_start[s]+
						  l*numsites+j][k], 
						num_smooth*numsamples, MEAN));
		      else
			fprintf(p, "%f ", 
				exp(find_statistics(x[series_state_start[s]+
						      l*numsites+j][k], 
						    num_smooth*numsamples, MEAN)));

		    bool done=false;
		    for(m=0;m<=meas_smooth[k].num_measurements && !done;m++)
		      if(m==meas_smooth[k].num_measurements ||
			 (meas_smooth[k].meanval[m]!=MISSING_VALUE &&
			  meas_smooth[k].site[m]==(int)j &&
			  meas_smooth[k].serie_num[m]==(int)s))
			{
			  if(!ser[s].pr->is_log)
			    fprintf(p,"%f ",
				    (m<meas_smooth[k].num_measurements) ? 
				    meas_smooth[k].meanval[m] : MISSING_VALUE);
			  else
			    fprintf(p,"%f ",
				    (m<meas_smooth[k].num_measurements) ? 
				    exp(meas_smooth[k].meanval[m]) : MISSING_VALUE);
			  done=true;
			}
		    if(ser[s].comp_meas)
		      fprintf(p," -10000000");
		    fprintf(p,"\n");
		  }
		
		if(ser[s].comp_meas)
		  for(k=0;k<ser[s].comp_len;k++)
		    if(ser[s].comp_meas[k].site==(int)j)
		      {
			if(ser[s].meas[0].dt!=NoHydDateTime)
			  fprintf(p,"%s ",  ser[s].comp_meas[k].dt.syCh(5));
			else
			  fprintf(p,"%f ", ser[s].comp_meas[k].tm);
			
			for(l=0;l<=ser[s].numlayers;l++)
			  fprintf(p, "-10000000 ");
			
			if(!ser[s].pr->is_log)
			  fprintf(p,"%f\n",  
				  ser[s].comp_meas[k].meanval);
			else
			  fprintf(p,"%f\n", 
				  ser[s].comp_meas[k].meanval!=MISSING_VALUE ?
				  exp(ser[s].comp_meas[k].meanval) : 
				  MISSING_VALUE);
		      }
		
		if(!*filestart)
		  pclose(p);
		else
		  fclose(p);
	      }
	  }
    }
  
  
  if(plot_repar)
    {
      char **repar_name=new char*[numpar];
      for(i=0;i<numpar;i++)
	{
	  repar_name[i]=new char[300];
	  snprintf(repar_name[i],250, "Re-parametrized %s",
		  par_name_orig[i]);
	}
      
      // show each parameter:
      for(i=0;i<numpar;i++)
	show_parameter(parsample_repar[i], numsamples, repar_name[i], silent, 
		       filestart);
      
      // show mean, median and 95% credibility interval for
      // each parameter:
      for(j=0;j<(int)numpar;j++)
	printf("%s: mean=%f median=%f    95%% cred=(%f-%f)\n", 
	       repar_name[j], 
	       find_statistics(parsample_repar[j], numsamples, MEAN),
	       find_statistics(parsample_repar[j], numsamples, MEDIAN),
	       find_statistics(parsample_repar[j], numsamples, PERCENTILE_2_5),
	       find_statistics(parsample_repar[j], numsamples, PERCENTILE_97_5));
      
      // If not silent, show also scatter diagrams for
      // each combination of parameters:
      if(!silent)
	{
	  for(i=0;i<numpar;i++)
	    for(j=i+1;j<(int)numpar;j++)
	      show_scatter(parsample_repar[i], parsample_repar[j], numsamples, 
			   repar_name[i], repar_name[j], filestart);
	}
      
      doubledelete(repar_name,numpar);
    }
  
  
  if(plot_sd || report_sd)
    {
      for(s=0;s<num_series;s++)
	{
	  for(l=0;l<ser[s].numlayers;l++)
	    {
	      bool all_global=false;
	      for(j=0;j<(int)numsites && !all_global;j++)
		{
		  int  index_diffusion=-1, index_pull=-1;
		  bool found_diffusion=false, found_pull=false;
		  bool global_diffusion=false, global_pull=false;
		  
		  for(i=0;i<numpar && (!found_diffusion || !found_pull);i++)
		    {
		      if(par_series[i]==(int)s && par_layer[i]==(int)(l+1) &&
			 (par_region[i]==(int)j || par_region[i]<0))
			{
			  if(par_type[i]==SIGMA)
			    {
			      found_diffusion=true;
			      index_diffusion=i;
			      if(par_region[i]<0)
				global_diffusion=true;
			    }
			  else if(par_type[i]==DT)
			    {
			      found_pull=true;
			      index_pull=i;
			      if(par_region[i]<0)
				global_pull=true;
			    }
			}
		    }
		  
		  if(global_diffusion && global_pull)
		    all_global=true;
		  
		  if(found_diffusion && found_pull)
		    {
		      double *sd=new double[numsamples];
		      char s_name[1000];
		      
		      if(all_global)
			snprintf(s_name,999, "sd_s%d_l%d", s+1,l+1);
		      else
			snprintf(s_name,999, "sd_s%d_l%d_r%d", s+1,l+1,j);

		      for(i=0;(int)i<numsamples;i++)
			sd[i]=parsample[index_diffusion][i]*
			  sqrt(parsample[index_pull][i]/2.0);
		      
		      printf("%s: mean=%f median=%f    95%% cred=(%f-%f)\n", 
			     s_name,
			     find_statistics(sd, numsamples, MEAN),
			     find_statistics(sd, numsamples, MEDIAN),
			     find_statistics(sd, numsamples, PERCENTILE_2_5),
			     find_statistics(sd, numsamples, PERCENTILE_97_5));
		      
		      show_parameter(sd, numsamples, s_name, silent, filestart);
		      
		      delete [] sd;
		    }
		  else
		    {
		      if(!found_pull)
			cout << "Couldn't find pull for series " << s << 
			  " layer " << l+1 << " site " << j << endl;
		      else
			cout << "Couldn't find diffusion for series " << s << 
			  " layer " << l+1 << " site " << j << "!" << endl;
		    }
		}
	    }
	}
    }

  // show each parameter:
  for(i=0;i<numpar;i++)
    show_parameter(parsample[i], numsamples, par_name[i], silent, filestart);
  
  // If not silent, show also scatter diagrams for
  // each combination of parameters:
  if(!silent)
    {
      for(i=0;i<numpar;i++)
	for(j=i+1;j<(int)numpar;j++)
	  show_scatter(parsample[i], parsample[j], numsamples, 
		       par_name[i], par_name[j], filestart);
    }
  
  for(j=0;j<(int)numpar;j++)
    if(!strncasecmp(par_name[j],"beta",4))
      {
	int belowzero=0;
	for(i=0;(int)i<numsamples;i++)
	  if(parsample[j][i]<0.0)
	    belowzero++;
	cout << par_name[j] << " " << 
	  double(belowzero)/double(numsamples)*100.0 
	     << " % below zero" << endl;
      }
  
  reset_global_variables();
  delete [] use_site;
  doubledelete(parsample,numpar);
  doubledelete(parsample_repar,numpar);
  delete [] pars;
}  // END 

#endif // MAIN
