#include <stats.h>
#include <stdio.h>
#include <qfloat.h>
/*
10,5 0.01 = 0.37256530511457266
10,5 0.5 = 0.67424884471378727
10,5 0.99 = 0.89807140427814909
*/
int main(void)
{
        long double b = beta_distribution_inv(10.0 ,5.0 ,0.01 );
        printf("10,5 0.01 = %.19Lg\n",b);
        b = beta_distribution_inv(10.0 ,5.0 ,0.5);
        printf("10,5 0.5 = %.19Lg\n",b);
        b = beta_distribution_inv(10.0 ,5.0 ,0.99 );
        printf("10,5 0.99 = %.19Lg\n",b);

	b = beta_distribution(10.0,5.0,0.89807140427814909);
	printf("1/10,5,0.89807140427814909=%.19Lg\n",b);


	qfloat q = beta_distribution_invq(10.0q,5.0q,0.01q);
	printf("10,5 0.01 = %.79qg\n",q);
	q = beta_distribution_invq(10q ,5q ,0.5q);
	printf("10,5 0.5 = %.106qg\n",q);
	q = beta_distribution_invq(10.0q,5.0q ,0.99q);
	printf("10,5 0.99 = %.79qg\n",q);


	q = beta_distribution_invq(10q ,5q ,0.5q);
	q = beta_incompleteq(10q,5q,q);
	printf("\n*\n10,5 0.5 = %.79qg\n",q);

        q = 0.674248844713787273126479873997280054527971923977407208725236151620726650213422655301368619186127510663745Q;
        q = beta_incompleteq(10q,5q,q);
        printf("\n*\n10,5 0.5 = %.79qg\n",q);

	qfloat p;
	p = beta_distribution_invq(10q,5q,0.5q);
	q = beta_distribution_invq(5q,10q,0.5q);
	printf("10,5,0.5=\n%.105qg\n5,10,0.5=\n%.105qg\np+q=\n%.105qg\n",p,q,p+q);

}
