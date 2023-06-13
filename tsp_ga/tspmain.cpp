/*
The GAlib-based genetic algorithm code for the Travelling Salesman Problem (TSP) Berlin52.
���ߣ�wying
��λ����������ѧ���ѧԺ
*/

#include <math.h>
#include <ctime>
#include "ga/GASStateGA.h"
#include "ga/GASimpleGA.h"
#include "ga/GA1DArrayGenome.h"
#include "ga/garandom.h"
#include "ga/std_stream.h"

#define cout STD_COUT
#define cerr STD_CERR
#define endl STD_ENDL
#define ostream STD_OSTREAM
#define ifstream STD_IFSTREAM

// Set this up for your favorite TSP.  The sample one is a contrived problem
// with the towns laid out in a grid (so it is easy to figure out what the 
// shortest distance is, and there are many different paths with the same
// shortest path).  File format is that used by the TSPLIB problems.  You can 
// grab more problems from TSPLIB.
// 
#define MAX_TOWNS 100
#define TSP_FILE "berlin52.txt"

int x[MAX_TOWNS],y[MAX_TOWNS];//ÿ�����е�x�����y����
int DISTANCE[MAX_TOWNS][MAX_TOWNS];//ÿ��������֮������гɱ����ǶԳƵ�


float TSPObjective(GAGenome&);//����Ⱦɫ��������ܷ��õ�Ŀ�꺯��
void  TSPInitializer(GAGenome&);//TSP�����Ⱦɫ���ʼ������
int   TSPMutator(GAGenome&, float);//���TSP�����Ⱦɫ���������
int   TSPCrossover(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);//���TSP�����Ⱦɫ�彻������

void writeTSPPath(ostream & os, GAGenome& g);//��ָ��Ⱦɫ�������·�������ָ���ļ�

int main() {
  cout << "The GAlib program for the Travelling Salesman Problem (TSP) Berlin52.\n" << endl;

  //��Berlin52.txt�ļ���������������
  double CityID;
  ifstream in(TSP_FILE); 
  if(!in) {
    cerr << "could not read data file " << TSP_FILE << "\n";
    exit(1);
  }
  int ntowns=0;
  do {
    in >> CityID;
    in >> x[ntowns];
    in >> y[ntowns];
    ntowns++;
  } while(!in.eof() && ntowns < MAX_TOWNS);
  in.close();
  if(ntowns >= MAX_TOWNS) {
    cerr << "data file contains more towns than allowed for in the fixed\n";
    cerr << "arrays.  Recompile the program with larger arrays or try a\n";
    cerr << "smaller problem.\n";
    exit(1);
  }

  //���������������м�����гɱ�
  double dx,dy;
  for(int i=0;i<ntowns;i++) {
    for(int j=i; j<ntowns;j++) {
      dx=x[i]-x[j]; dy=y[i]-y[j];
      DISTANCE[j][i]=DISTANCE[i][j]=floor(0.5 + sqrt(dx*dx+dy*dy) );//ע��ȡ��������֮�������ֵ
    }
  }

  //����TSP����ı��뷽��Ϊһά���������飬��̶�����Ϊ���и���
  GA1DArrayGenome<int> genome(ntowns);

  genome.evaluator(::TSPObjective);//ΪȾɫ��ָ������Ŀ��ֵ�ĺ���
  genome.initializer(::TSPInitializer);//ΪȾɫ��ָ���Զ���ĳ�ʼ������
  genome.crossover(::TSPCrossover);//ΪȾɫ��ָ���Զ���Ľ�������
  genome.mutator(::TSPMutator);//ΪȾɫ��ָ���Զ���ı�������

  GASteadyStateGA ga(genome); ga.nReplacement(2); ga.nGenerations(500000);//ѡ����̬�Ŵ��㷨����TSP������⣬ָ����Ⱦɫ����뷽ʽ��ÿһ��Ҫ�滻�ĸ�����=2���ܵ����д���500000����ô�������ܸ�����=2*500000=1000000
  //GASimpleGA ga(genome); ga.elitist(gaTrue); ga.nGenerations(5000);//ѡ�ü��Ŵ��㷨����TSP������⣬���þ�Ӣ�������ԣ�ָ����Ⱦɫ����뷽ʽ���ܵ����д���10000����ô�������ܸ�����=200����Ⱥ��С��*5000=1000000
  ga.minimize();//Ϊ�Ŵ��㷨ָ���Ż�Ŀ���ǽ�Ŀ�꺯��ֵ��С��
  ga.populationSize(200);//Ϊ�Ŵ��㷨ָ����Ⱥ��СΪ200
  ga.pMutation(0.02);//Ϊ�Ŵ��㷨ָ���������
  ga.pCrossover(0.8);//Ϊ�Ŵ��㷨ָ���������

  cout << "initializing..."<<"\n"; cout.flush();
  unsigned int seed = clock();
  ga.initialize(seed);//ʹ�ô�ʱ�ӵõ���������ӳ�ʼ���Ŵ��㷨

  cout << "evolving..."<<"\n"; cout.flush();
  std::fstream fgacurve;
  fgacurve.open("tspgacurve.txt",std::ios::out);

  //�Ŵ��㷨��ʼ����������ֱ���ﵽָ���Ĵ���
  while(!ga.done()) {
    ga.step();//����һ��
    if(ga.generation() % (ga.nGenerations()/100) == 0) 
	{//����������ȡ100�������㣬��¼���������е�����Ŀ��ֵ������Ϣ���ļ�
		int bestscore = ga.statistics().bestIndividual().score();
		cout << ga.generation() << "    " << bestscore << "\n"; cout.flush();
		fgacurve << ga.generation() << "    " << bestscore << "\n";
    }
  }
  fgacurve.close();

  //�Ŵ��㷨������ֹ������ҵ�����������·�ߵ��ļ�
  genome = ga.statistics().bestIndividual();
  //cout << "\n" << "the shortest path found is "  << "\n";
  //writeTSPPath(cout, genome);
  std::fstream fbestpath;
  fbestpath.open("tsppath.txt",std::ios::out);
  writeTSPPath(fbestpath, genome);
  fbestpath.close();
  cout << "the distance of the shortest path found: " << genome.score() << "\n";

  return 0;
}


// Here are the genome operators that we want to use for this problem.
//����Ⱦɫ��������ܷ��õ�Ŀ�꺯��
float TSPObjective(GAGenome& g) {
  GA1DArrayGenome<int> & genome = (GA1DArrayGenome<int> &)g;
  int genomelength = genome.size();//genome.size()��ȡȾɫ��ĳ���
  float dist = 0;
  int xx;
  int yy;

    for(int i=0; i<genomelength; i++) {
      xx = genome.gene(i);
      yy = genome.gene( (i+1)%genomelength );
      dist += DISTANCE[xx-1][yy-1];
    }

  return dist;
}

//TSP�����Ⱦɫ���ʼ������
void TSPInitializer(GAGenome& g) {
  GA1DArrayGenome<int> &genome=(GA1DArrayGenome<int> &)g;

  int genomelength = genome.size();
  int i,town;
  static bool visit[MAX_TOWNS];

  memset(visit, false, MAX_TOWNS*sizeof(bool));
  town=GARandomInt(1,genomelength);//GARandomInt(1,genomelength)����1��genomelength֮���һ�������������
  visit[town-1]=true;
  genome.gene(0, town);//genome.gene(0, town)���ø�Ⱦɫ���0������λ�ϵĻ���ֵΪtown
 
  for( i=1; i<genomelength; i++) {
    do {
      town=GARandomInt(1,genomelength);
    } while (visit[town-1]);
    visit[town-1]=true;
    genome.gene(i, town);
  }	
}

//���TSP�����Ⱦɫ��������ӣ�pmutΪ�������
int TSPMutator(GAGenome& g, float pmut) {
  GA1DArrayGenome<int> &genome=(GA1DArrayGenome<int> &)g;
  int i;

  int genomelength = genome.size();
  float nmutator = pmut*genomelength;//Ҫ�ı�ıߵ�����

  int imutator=0;
  while( imutator<nmutator){
    if (GARandomFloat()<0.5) {//GARandomFloat()����0��1֮���һ���������������
	  //��0.5����ʹ���໥��������
	  int swapIndex1 = GARandomInt(0,genomelength-1);
	  int swapIndex2 = GARandomInt(0,genomelength-1);
	  int tmp;
	  tmp = genome.gene(swapIndex2);
	  genome.gene(swapIndex2, genome.gene(swapIndex1) );
	  genome.gene(swapIndex1, tmp );// swap only one time
	  imutator+=4;
    }else
	  {
	  //��0.5����ʹ�÷�ת����
	  int inversion_start, inversion_end, tmp;
	  inversion_start = GARandomInt(0,genomelength-1);
	  inversion_end = GARandomInt(0,genomelength-1);
	  if(inversion_start > inversion_end)
	  {
		  tmp = inversion_start;
		  inversion_start = inversion_end;
		  inversion_end = tmp;
	  }
	   
	  for(i=inversion_start;inversion_start<inversion_end; inversion_start++, inversion_end-- )
	  {
		  tmp = genome.gene(inversion_start);
		  genome.gene(inversion_start, genome.gene(inversion_end) );
		  genome.gene(inversion_end, tmp ); 
	  }  
	  imutator+=2;
    }
  }

  return (1);
}


//���TSP�����Ⱦɫ�彻������
int TSPCrossover(const GAGenome& g1, const GAGenome& g2, GAGenome* c1, GAGenome* c2) {
  GA1DArrayGenome<int> parent1=(GA1DArrayGenome<int> &)g1;
  GA1DArrayGenome<int> parent2=(GA1DArrayGenome<int> &)g2;

  int genomelength = parent1.size();

  int nc=0;

  GA1DArrayGenome<int> &child1=(GA1DArrayGenome<int> &)*c1;
  GA1DArrayGenome<int> &child2=(GA1DArrayGenome<int> &)*c2;

  if(c1)  {child1 = parent2; nc++;}
  if(c2)  {child2 = parent1; nc++;}

  /////

  //cout <<"ԭ��1��"<< parent1 << endl <<"ԭ��2��"<< parent2 << endl; 
  //cout << genomelength << endl; ��鳤��
  /*int swapIndex1 = GARandomInt(0, genomelength - 1);
  int swapIndex2 = GARandomInt(0, genomelength - 1);  ģ��������������ķ���*/
  int swappoint1 = 0;
  int swappoint2 = 0;
  int temp = 0;
  do {
	  swappoint1= GARandomInt(1, genomelength - 2); //��һ��������
	  swappoint2 = GARandomInt(1, genomelength - 2); //�ڶ���������
  } while (swappoint1 == swappoint2);
  if (swappoint1 > swappoint2) {
	  temp = swappoint1;
	  swappoint1 = swappoint2;
	  swappoint2 = temp;
  }//ȷ����1�ڵ�2ǰ
  //cout << swappoint1 << "      " << swappoint2 << endl;
  int check1 = 0;//������
  int check2 = 0;
  int temp1 = 0;//�����Ƿ����
  int temp2 = 0;
  for (int i = 0; i <= genomelength - 1; i++) //������ɨ�裬��������뽻�����в�һ�µı��
  {
	  if (i == swappoint1) //
	  { 
		  i = swappoint2 + 1; //�����Ѿ������Ƭ��
	  }
	  while (check1 <= genomelength - 1)  //ѭ�����
	  {
		  temp1 = 0;
		  for (int k = swappoint1; k <= swappoint2; k++) 
		  {
			  if (parent1.gene(check1)==parent2.gene(k)) //�ͽ���Ƭ�ν��бȽϣ�һ�µ���������һ�µ����
			  {
				  temp1 = 1;
				  
				  break;
			  }
		  }
		  if (temp1==0)
		  {
			  child1.gene(i, parent1.gene(check1)); //���
			  check1 = check1 + 1;
			  break;

		  }
		  check1 = check1 + 1;
	  }

  }

  for (int i = 0; i <= genomelength - 1; i++)//ͬ��
  {
	  if (i == swappoint1)
	  {
		  i = swappoint2 + 1;
	  }
	  while (check2 <= genomelength - 1)
	  {
		  temp2 = 0;
		  for (int k = swappoint1; k <= swappoint2; k++)
		  {
			  if (parent2.gene(check2) == parent1.gene(k))
			  {
				  temp2 = 1;
				  
				  break;
			  }
		  }
		  if (temp2 == 0)
		  {
			  child2.gene(i, parent2.gene(check2));
			  check2 = check2 + 1;
			  break;

		  }
		  check2 = check2 + 1;
	  }

  }
  //cout << "����1��" << child1 << endl;
  //cout<< "����2��" << child2 << endl;
	//�˴���Ӵ���ʵ���Լ��Ľ�������
  /////
  return nc;
}

//��ָ��Ⱦɫ�������·�������ָ���ļ�
void writeTSPPath(ostream & os, GAGenome& g) 
{
	GA1DArrayGenome<int> & genome = (GA1DArrayGenome<int> &)g;
	int genomelength = genome.size();
	for(int i=0; i<genomelength; i++)
	{
		int xx = genome.gene(i);
		os << xx <<"    " <<x[xx-1] << "      "<<y[xx-1] << "\n";
	}
}
