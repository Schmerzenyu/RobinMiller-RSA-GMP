#include <stdio.h>
#include "iostream"
#include "gmp.h"
#include "ctime"
#pragma comment(linker, "/NODEFAULTLIB:MSVCRTD.LIB")
#pragma comment(lib, "gmp.lib")
#pragma comment(lib, "gmpDebug.lib")
using namespace std;

const int TestCase = 5;
#define PRIME_PROBABILITY 20
#define PRIME_BITLENGTH 1024

void CreateBigPrime(mpz_t mpzPrime,int bits)//,int Probability)
{
	int i,last_rand=0;
	char *char_rand = new char [bits+1];
	char_rand[0] = '1';
	char_rand[bits] = '\0';
	mpz_init(mpzPrime);

	do
	{
		for(i=1;i<bits;i++)
		{
			if(rand()==last_rand)
				cout<<"SRAND ERROR! --short time! same rand number!"<<endl;
			last_rand = rand();
			char_rand[i] = '0'+(0x01&last_rand);
		}
			
		mpz_set_str(mpzPrime,char_rand,2);
		mpz_nextprime(mpzPrime,mpzPrime);
	}while(0==mpz_probab_prime_p(mpzPrime,PRIME_PROBABILITY));
	
	return;
}


//a ^ b % n
void BigIntergerMod(mpz_t a, mpz_t b, mpz_t n, mpz_t &res)
{
	mpz_t temp;
	mpz_t k;
	mpz_t bb;
	mpz_init(temp);
	mpz_init(k);
	mpz_init(bb);
	mpz_set(bb,b);//bb = b;
	//res = 1;
	mpz_set_ui(res,1);
	mpz_mod(temp,a,n);//temp = a % n;
	while (mpz_cmp_ui(bb,0))//bb > 0 ?
	{
		mpz_mod_ui(k, bb, 2);// k = bb % 2;
		if (mpz_cmp_ui(k, 0))//k == 1?
		{
		//	res = res * temp % n;
			mpz_mul(res, res, temp);
			mpz_mod(res, res, n);
		}
		//bb = bb >> 1;
		mpz_fdiv_q_ui(bb, bb, 2);
		//temp = temp * temp % n;
		mpz_mul(temp, temp, temp);
		mpz_mod(temp,temp, n);
		
	}
	mpz_clear(temp);
	mpz_clear(k);
	mpz_clear(bb);
}


//判断n是否为素数，若为合数返回true，若可能为素数返回false
bool TestPrime(mpz_t n)
{
	mpz_t a, k, q, j, judge, n_sub, res;
	mpz_init(a);
	mpz_init(k);
	mpz_init(q);
	mpz_init(j);
	mpz_init(judge);
	mpz_init(n_sub);
	mpz_init(res);
	gmp_randstate_t state;
	gmp_randinit_default(state);
	mpz_sub_ui(n_sub, n, 1);//n_sub = n - 1;
	mpz_set(q, n_sub);//q = n - 1;
	mpz_mod_ui(judge, q, 2);// judge = q % 2;
	// n-1 = 2^k*q
	while (!mpz_cmp_ui(judge,0))//judge = 0
	{
		mpz_fdiv_q_ui(q, q, 2);//q = q / 2;
		mpz_add_ui(k, k, 1);// k++;
		mpz_mod_ui(judge, q, 2);// judge = q % 2;
	}
//	gmp_printf("k=%Zd\n",k);
//	gmp_printf("q=%Zd\n",q);
	mpz_urandomm(a,state,n);//产生随机数1 < a <= n-1;
	mpz_sub_ui(a, a, 1);//a = a - 1;
	BigIntergerMod(a, q, n, res);//res = a ^ q % n
	if (!mpz_cmp_ui(res, 1)||!mpz_cmp(res, n_sub))//res == 1||res == n - 1
		return  false;
	else
	{
		mpz_add_ui(j,j,1);//j++;
		for ( ;mpz_cmp(j, k); mpz_add_ui(j, j, 1))//j = 1..k-1
		{
			mpz_mul(res, res, res);//res = res ^ 2;
			mpz_mod(res, res, n);//res = res % n;
			if ( !mpz_cmp_ui(res, 1))//res == 1
				return true;
			if ( !mpz_cmp(res, n_sub) )// res == n-1
				return false;
		}
	}
	mpz_clear(a);
	mpz_clear(k);
	mpz_clear(q);
	mpz_clear(j);
	mpz_clear(judge);
	mpz_clear(n_sub);
	mpz_clear(res);
	return true;
}

//Miller Rabin素数测试，是素数返回true，否则返回false
bool MRTest(mpz_t n)
{
	mpz_t temp;
	mpz_init(temp);
	mpz_mod_ui(temp, n, 2);
	if (!mpz_cmp_ui(temp, 0))//temp == 0
	{
		mpz_clear(temp);
		return false;
	}
	mpz_clear(temp);
	for ( int i = 0; i < TestCase; ++i )
	{
		if (TestPrime(n))
			return false;
	}
	return true;
}



int main()
{   
	mpz_t p,q,e,d,p_sub,q_sub,n,z,temp;
	CreateBigPrime(p, PRIME_BITLENGTH);
	CreateBigPrime(q, PRIME_BITLENGTH);
	gmp_printf("生成的素数\np= %Zd\nq= %Zd\n",p,q);
//  n =pq z=(p-1)(q-1)
	mpz_init(n);
	mpz_init(p_sub);
	mpz_init(q_sub);
	mpz_init(z);
	mpz_mul(n, p, q);
	mpz_sub_ui(p_sub, p, 1);
	mpz_sub_ui(q_sub, q, 1);
	mpz_mul(z, p_sub, q_sub);
	gmp_printf("乘积\nn= %Zd\nz= %Zd\n",n,z);
//  选取e
	mpz_init(temp);
	mpz_init(e);
	int cmp_result = 1;
	for(int i = 2; cmp_result != 0; i=i+1 )
	{
		mpz_gcd_ui(temp, z, i);
		cmp_result = mpz_cmp_ui (temp, 1) ;
		if(cmp_result ==0)
		{
		mpz_init_set_ui(e, i);
		}	
	}
//  得到逆元d
	mpz_init(d);
	mpz_invert (d, e, z);
	gmp_printf("随机互素的e为\ne= %Zd\n逆元\nd= %Zd\n",e,d);
//	公钥为(n,e)
	gmp_printf("公钥为\nn= %Zd\nd= %Zd\n",e,d);

//	私钥为(p,q,d)
	gmp_printf("私钥为\np= %Zd\nq= %Zd\nd= %Zd\n",p,q,d);
	
	//mpz_init(key);
	mpz_t content;
	char *keystr=new char[100];
	printf("请输入你要加密的明文（数字，100位以内）");
	gets(keystr);

	mpz_init_set_str(content,keystr, 10);
    gmp_printf("明文为%Zd\n",content);
	//m ^ e % n
	mpz_t secret;
	mpz_init(secret);
	BigIntergerMod(content, e, n, secret);
	gmp_printf("密文为%Zd\n",secret);
	//c ^ d % n
	BigIntergerMod(secret, d, n, content);
	gmp_printf("解密得到的明文为%Zd\n", content);
	mpz_clear(p);
	mpz_clear(q);
	mpz_clear(p_sub);
	mpz_clear(q_sub);
	mpz_clear(e);
	mpz_clear(n);
	mpz_clear(z);
	mpz_clear(d);
	mpz_clear(content);
	mpz_clear(secret);
	system("pause"); 
	return 0;
}
