#ifndef _RBF_GENOME_HEADER_
#define _RBF_GENOME_HEADER_

#include <iostream>
#include <ga/ga.h>

#define INSTANTIATE_REAL_GENOME
#include <ga/GARealGenome.h>

class RBFGenome : public GAGenome {
public:
  GADefineIdentity("RBFGenome", 201);

  static void RBFInitializer(GAGenome&);
  static int RBFMutator(GAGenome&, float);
  static float RBFComparator(const GAGenome&, const GAGenome&);
  static int RBFCrossover(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);
public:
  RBFGenome(int,GAAlleleSet<float>&,GAGenome::Evaluator f=NULL, void* u=NULL);
  RBFGenome(const RBFGenome & orig);
  RBFGenome& operator=(const GAGenome& g);
  virtual ~RBFGenome();
  virtual GAGenome* clone(GAGenome::CloneMethod) const ;
  virtual void copy(const GAGenome & c);
  virtual int equal(const GAGenome& g) const;
  virtual int read(std::istream & is);
  virtual int write(std::ostream & os) const;

  GA1DBinaryStringGenome & binstr() const {return *str;}
  GARealGenome & real() const {return *r;}

protected:
  GA1DBinaryStringGenome *str;
  GARealGenome *r;
};


//// Member functions for the RBF genome object
RBFGenome::
RBFGenome(int length, GAAlleleSet<float>& alleleset, GAGenome::Evaluator f, void* u) :
		GAGenome(  RBFInitializer,
			 RBFMutator,
			 RBFComparator) {
		  evaluator(f); userData(u); crossover(RBFCrossover);

		  str = new GA1DBinaryStringGenome(length, f, u);
		  r = new GARealGenome(length, alleleset, f, u);
		}

RBFGenome::RBFGenome(const RBFGenome & orig) {
  str = new GA1DBinaryStringGenome(orig.binstr());
  r = new GARealGenome(orig.real());
  copy(orig);
}

RBFGenome&
RBFGenome::operator=(const GAGenome& g) { copy(g); return *this; }

RBFGenome::~RBFGenome() { delete str; delete r; }

GAGenome*
RBFGenome::clone(GAGenome::CloneMethod) const {
  return new RBFGenome(*this);
}

void
RBFGenome::copy(const GAGenome & c){
  if(&c != this && sameClass(c)){
    GAGenome::copy(c);
    RBFGenome & bc = (RBFGenome &)c;
    str->copy(*(bc.str));
    r->copy(*(bc.r));
  }
}

int
RBFGenome::equal(const GAGenome& g) const {
  RBFGenome& genome = (RBFGenome&)g;
  return ((str->equal(*genome.str)) && (r->equal(*genome.r)));
}

int
RBFGenome::read(std::istream & is) {
  is >> *str >> *r;
  return is.fail() ? 1 : 0;
}

int
RBFGenome::write(std::ostream & os) const {
  os << *str <<"\n" << *r << "\n";
  return os.fail() ? 1 : 0;
}



// These are the default initialization, mutation, and comparator operators for
// this genome class.  They are defined as static functions of the RBF
// genome class and they're defaults for the class.  But they can be overridden
// on any instance of the genome.

// The initializer just calls the initializer for each of the genomes that are
// in the RBF genome.

// I would have used simply 'Initializer', 'Mutator', etc rather than
// 'RBFInitializer' but old versions of g++ are brain-dead and don't
// get the encapsulation properly.
void
RBFGenome::RBFInitializer(GAGenome & c) {
  RBFGenome & child = (RBFGenome &)c;
  child.binstr().initialize();
  child.real().initialize();
  child._evaluated = gaFalse;
}


// The mutator just calls the mutator for each of the component genomes.
int
RBFGenome::RBFMutator(GAGenome & c, float pmut) {
  RBFGenome & child = (RBFGenome &)c;
  int nmut = child.binstr().mutate(pmut) + child.real().mutate(pmut);
  if(nmut) child._evaluated = gaFalse;
  return nmut;
}

// The comparator just calls the comparators for each of the component genomes,
// then averages the score.
float
RBFGenome::RBFComparator(const GAGenome& a, const GAGenome& b) {
  RBFGenome& sis = (RBFGenome &)a;
  RBFGenome& bro = (RBFGenome &)b;
  return 0.5 * (sis.binstr().compare(bro) + sis.real().compare(bro));
}

// The crossover operator invokes the crossover for each of the genomes in the
// RBF genome.  We use sexual crossover only, and we do not test to see
// if no crossover has been assigned.
int
RBFGenome::
RBFCrossover(const GAGenome& a, const GAGenome& b,
		   GAGenome* c, GAGenome* d){
  RBFGenome& mom = (RBFGenome&)a;
  RBFGenome& dad = (RBFGenome&)b;
  int n=0;

  GAGenome::SexualCrossover strcross = mom.str->sexual();
  GAGenome::SexualCrossover rcross = mom.r->sexual();

  if(c && d){
    RBFGenome& sis = (RBFGenome&)*c;
    RBFGenome& bro = (RBFGenome&)*d;
    (*strcross)(mom.binstr(), dad.binstr(), &sis.binstr(), &bro.binstr());
    (*rcross)(mom.real(),dad.real(), &sis.real(), &bro.real());
    sis._evaluated = gaFalse;
    bro._evaluated = gaFalse;
    n = 2;
  }
  else if(c){
    RBFGenome& sis = (RBFGenome&)*c;
    (*strcross)(mom.binstr(), dad.binstr(), &sis.binstr(), 0);
    (*rcross)(mom.real(), dad.real(), &sis.real(), 0);
    sis._evaluated = gaFalse;
    n = 1;
  }
  else if(d){
    RBFGenome& bro = (RBFGenome&)*d;
    (*strcross)(mom.binstr(), dad.binstr(), 0, &bro.binstr());
    (*rcross)(mom.real(), dad.real(), 0, &bro.real());
    bro._evaluated = gaFalse;
    n = 1;
  }

  return n;
}

#endif
