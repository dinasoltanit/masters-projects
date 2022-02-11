#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void)
{
    /*
    Dear user,
    you can take your data by asking the user in method 1
    or you can simply edit the following lines using method2.
    Each time that you run the code, you will be asked these
    general questions:
    1. Please specify the method you want to work with.
    2. Please indicate the percentage of the grades' contribution.
    All the best
    First, you need to choose your method.
    Enter 1 if you want method 1, and 2 if you want method 2.
    Method 1:
    Enter 1 if you want to choose this method.
    In this method, you'll go through a loop. In each loop,
    you are asked to insert the name of the student, following
    with the grades from each of the sections.
    Method 2:
    Enter 2 if you want to choose this method. In this method,
    you need to simply modify the data below, based on your own
    data.
    */
    int Method = 2;
    float percents[7] = {0.1, 0.06, 0.08, 0.08, 0.08, 0.2, 0.4};
    //char input[50];
    char M1_name[10][2][20];
    float M1_HW[10];
    float M1_P1[10];
    float M1_P2[10];
    float M1_P3[10];
    float M1_P4[10];
    float M1_mid[10];
    float M1_fin[10];
    int M1_Absences[10];
    int M1_Delays[10];
    float M1_Grade[10];

    char M2_name[10][2][8] = {{"Mina    ", "Seyfan  "}, {"Ali     ", "Hosseini"}, {"Setare  ", "Iraji   "}, {"Parsa   ", "Seyedi  "}, {"Fateme  ", "Keshtkar"}, {"Mostafa ", "Kia     "}, {"Maryam  ", "Seif    "}, {"Bardia  ", "Naji    "}, {"Soheila ", "Sarlak  "}, {"Majid   ", "Hadi    "}};
    //char M2_name[10][20] = {"Mina Seyfan", "Ali Hosseini", "Setared Iraji", "Parsa Seyedi", "Fateme Keshtkar", "Mostafa Kia", "Maryam Seif", "Bardia Naji", "Soheila Sarlak", "Majid Hadi"};
    float M2_HW[10] = {91, 91.5, 92, 92.5, 93, 93.5, 94, 94, 94, 94};
    float M2_P1[10] = {91, 91.5, 92, 92.5, 93, 93.5, 94, 94, 94, 94};
    float M2_P2[10] = {91, 91.5, 92, 92.5, 93, 93.5, 94, 94, 94, 94};
    float M2_P3[10] = {91, 91.5, 92, 92.5, 93, 93.5, 94, 94, 94, 94};
    float M2_P4[10] = {91, 91.5, 92, 92.5, 93, 93.5, 94, 94, 94, 94};
    float M2_mid[10] = {91, 91.5, 92, 92.5, 93, 93.5, 94, 94, 94, 94};
    float M2_fin[10] = {91, 91.5, 92, 92.5, 93, 93.5, 94, 94, 94, 94};
    int M2_Absences[10] = {1, 2, 3, 5, 4, 0, 1, 2, 3, 4};
    int M2_Delays[10] = {1, 2, 3, 5, 4, 0, 1, 2, 3, 4};
    float M2_Grade[10];
    int i, j, k;

    printf("Please enter 1 to choose method1, and 2 to choose method2: ");
    scanf("%d", &Method);
    printf("You entered %d.\n", Method);

    switch (Method)
    {
    case 1:
        /* ******************************Case 1********************************** */
        printf("Enter the percent for the homeworks; the default is 0.1: ");
        scanf("%f", &percents[0]);
        printf("You entered %f.\n", percents[0]);
        printf("Enter the percent for the project number 1; the default is 0.06: ");
        scanf("%f", &percents[1]);
        printf("Enter the percent for the project number 2; the default is 0.08: ");
        scanf("%f", &percents[2]);
        printf("Enter the percent for the project number 3; the default is 0.08: ");
        scanf("%f", &percents[3]);
        printf("Enter the percent for the project number 4; the default is 0.08: ");
        scanf("%f", &percents[4]);
        printf("Enter the percent for the midterm exam; the default is 0.2: ");
        scanf("%f", &percents[5]);
        printf("Enter the percent for the final exam; the default is 0.4: ");
        scanf("%f", &percents[6]);
        for (i = 0; i < 10; i++)
        {
            printf("You are in loop number %d.\n", i);

            printf("Enter the given name of the student number %d: ", i);
            scanf("%s", &M1_name[i][0][20]);
            printf("Enter the last name of the student number %d: ", i);
            scanf("%s", &M1_name[i][1][20]);

            /*gets(input);
            char *token = strtok(input, " ");
            strcpy(M1_name[i][1][20], token);
            //printf("%s\n",command);
            token = strtok(NULL, " ");
            strcpy(M1_name[i][2][20], token);
            //printf("%s\n",num);*/

            printf("\nEnter the HW's mark out of 100: ");
            scanf("%f", &M1_HW[i]);
            printf("\nEnter the 1st Project's mark out of 100: ");
            scanf("%f", &M1_P1[i]);
            printf("\nEnter the 2nd Project's mark out of 100: ");
            scanf("%f", &M1_P2[i]);
            printf("\nEnter the 3rd Project's mark out of 100: ");
            scanf("%f", &M1_P3[i]);
            printf("\nEnter the 4th Project's mark out of 100: ");
            scanf("%f", &M1_P4[i]);
            printf("\nEnter the Midterm Exam's mark out of 100: ");
            scanf("%f", &M1_mid[i]);
            printf("\nEnter the Final Exam's mark out of 100: ");
            scanf("%f", &M1_fin[i]);
            printf("\nEnter the number of absences: ");
            scanf("%d", &M1_Absences[i]);
            printf("\nEnter the number of delays: ");
            scanf("%d", &M1_Delays[i]);
            M1_Grade[i] = M1_HW[i] * percents[0] + M1_P1[i] * percents[1] + M1_P2[i] * percents[2] + M1_P3[i] * percents[3] + M1_P4[i] * percents[4] + M1_mid[i] * percents[5] + M1_fin[i] * percents[6] + (4 - M1_Absences[i] - (M1_Delays[i] / 3));
            //printf("grade = ");
            //printf("grade = %f\n", M1_Grade[i]);
        }
        printf("Name    Grade\n");
        for (j = 0; j < 10; j++)
        {
            printf("%s %s   %f \n", &M1_name[j][0][20], &M1_name[j][1][20], M1_Grade[j]);
        }
        break;
    case 2:
        /* ******************************Case 2********************************** */
        printf("Enter the percent for the homeworks; the default is 0.1: ");
        scanf("%f", &percents[0]);
        printf("You entered %f.\n", percents[0]);
        printf("Enter the percent for the project number 1; the default is 0.06: ");
        scanf("%f", &percents[1]);
        printf("Enter the percent for the project number 2; the default is 0.08: ");
        scanf("%f", &percents[2]);
        printf("Enter the percent for the project number 3; the default is 0.08: ");
        scanf("%f", &percents[3]);
        printf("Enter the percent for the project number 4; the default is 0.08: ");
        scanf("%f", &percents[4]);
        printf("Enter the percent for the midterm exam; the default is 0.2: ");
        scanf("%f", &percents[5]);
        printf("Enter the percent for the final exam; the default is 0.4: ");
        scanf("%f", &percents[6]);

        printf("You are using method 2.\n");
        printf("Name    Grade\n");
        for (k = 0; k < 10; ++k)
        {
            //printf("Given name: %s, Last name: %s \n", &M2_name[0][0][0], &M2_name[0][1][0]);
            //printf("Given name: %s, Last name: %s \n", &M2_name[k][0][0], &M2_name[k][1][0]);
            M2_Grade[k] = M2_HW[k] * percents[0] + M2_P1[k] * percents[1] + M2_P2[k] * percents[2] + M2_P3[k] * percents[3] + M2_P4[k] * percents[4] + M2_mid[k] * percents[5] + M2_fin[k] * percents[6] + (4 - M2_Absences[k] - (M2_Delays[k] / 3));
            //printf("Calculating ... loop %d \n", k);

            printf("%s %s   %f \n", &M2_name[k][0], &M2_name[k][1], M2_Grade[k]);
            //printf("%s   %f \n", &M2_name[k][20], M2_Grade[k]);
        }
        break;
    default:
        printf("You are using the default setting.\n");
        printf("Name    Grade\n");
        for (k = 0; k < 10; ++k)
        {
            //printf("Given name: %s, Last name: %s \n", &M2_name[0][0][0], &M2_name[0][1][0]);
            //printf("Given name: %s, Last name: %s \n", &M2_name[k][0][0], &M2_name[k][1][0]);
            M2_Grade[k] = M2_HW[k] * percents[0] + M2_P1[k] * percents[1] + M2_P2[k] * percents[2] + M2_P3[k] * percents[3] + M2_P4[k] * percents[4] + M2_mid[k] * percents[5] + M2_fin[k] * percents[6] + (4 - M2_Absences[k] - (M2_Delays[k] / 3));
            //printf("Calculating ... loop %d \n", k);

            printf("%s %s   %f \n", &M2_name[k][0], &M2_name[k][1], M2_Grade[k]);
            //printf("%s   %f \n", &M2_name[k][20], M2_Grade[k]);
        }
        break;
    }
    return 0;
}
