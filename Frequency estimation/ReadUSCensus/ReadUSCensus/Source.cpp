#include <stdio.h>
#include <stdlib.h>
#include <windows.h>

#define INI_PATH		".\\ReadUSCensus.ini"
#define MAX_CAT_NUM		68

int ReadCatValNum[MAX_CAT_NUM];

struct iniParameter{
	char uscen_file[129];
	char catid_file[129];
	int read_cat[MAX_CAT_NUM];
	int read_cat_num;
};

void Exit(int rc){
	printf("\nPress ENTER to exit.\n"); getchar();
	exit(rc);
}

void FileOpen(FILE **file, char *filename, char *mode){
	errno_t error;

	if (error = fopen_s(file, filename, mode) != 0){
		printf("Error: Canot open %s\n", filename);
		Exit(-1);
	}
}

iniParameter readIniParameter(){
	iniParameter ret;
	char s[1025];
	char *tok, *ctx;

	GetPrivateProfileString("PATH", "USCEN_FILE", "default", ret.uscen_file, sizeof(ret.uscen_file), INI_PATH);
	GetPrivateProfileString("PATH", "CATID_FILE", "default", ret.catid_file, sizeof(ret.catid_file), INI_PATH);

	GetPrivateProfileString("PARAMETER", "READ_CAT", "default", s, sizeof(s), INI_PATH);
	tok = strtok_s(s, ",", &ctx);
	ret.read_cat_num = 0;
	ret.read_cat[ret.read_cat_num++] = atoi(tok);
	while ((tok = strtok_s(NULL, ",", &ctx)) != NULL){
		ret.read_cat[ret.read_cat_num++] = atoi(tok);
	}

	return ret;
}

void ReadUSCensus(char *uscen_file, int *read_cat, int read_cat_num, char *catid_file){
	int *cat_tmp;
	int tot_cat_id;
	int user_num;
	char s[10001];
	int i, j;
	FILE *fp, *fp2;
	char *tok, *ctx;

	cat_tmp = (int *)malloc(sizeof(int) * read_cat_num);

	FileOpen(&fp, uscen_file, "r");
	FileOpen(&fp2, catid_file, "w");

	fgets(s, 10000, fp);

	user_num = 0;
	printf("Extracting Category ID (#:100000 users): ");
	while (fgets(s, 10000, fp) != NULL){
		if (user_num % 100000 == 0) printf("#");
		if ((tok = strtok_s(s, ",", &ctx)) == NULL){
			printf("Error: File format is incorrect.\n");
			Exit(-1);
		}

		for (j = 0; j < read_cat_num; j++) cat_tmp[j] = 0;

		for (i = 0; i < 68; i++){
			if ((tok = strtok_s(NULL, ",", &ctx)) == NULL){
				printf("Error: File format is incorrect.\n");
				Exit(-1);
			}

			for (j = 0; j < read_cat_num; j++){
				if (i == read_cat[j] - 1){
					cat_tmp[j] = atoi(tok);
					break;
				}
			}
		}

		tot_cat_id = cat_tmp[0];
		for (j = 1; j < read_cat_num; j++){
			tot_cat_id *= ReadCatValNum[read_cat[j] - 1];
			tot_cat_id += cat_tmp[j];
		}

		fprintf_s(fp2, "%d,", user_num + 1);
		fprintf_s(fp2, "%d,", tot_cat_id);
		for (j = 0; j < read_cat_num; j++) fprintf_s(fp2, "%d,", cat_tmp[j]);
		fprintf_s(fp2, "\n");

		user_num++;
	}
	printf("done.\n");

	fclose(fp);
	fclose(fp2);
	free(cat_tmp);
}

void main(){
	struct iniParameter ini = readIniParameter();

	ReadCatValNum[0]  =  8;		// 1:dAge
	ReadCatValNum[5]  = 10;		// 6:iClass
	ReadCatValNum[7]  =  3;		// 8:iDisabl1
	ReadCatValNum[8]  =  3;		// 9:iDisabl2
	ReadCatValNum[9]  =  5;		// 10:iEnglish
	ReadCatValNum[10] =  2;		// 11:iFeb55
	ReadCatValNum[11] = 14;		// 12:iFertil
	ReadCatValNum[16] =  5;		// 17:dIncome1
	ReadCatValNum[25] =  2;		// 26:iKorean
	ReadCatValNum[26] =  3;		// 27:iLang1
	ReadCatValNum[27] =  3;		// 28:iLooking
	ReadCatValNum[28] =  5;		// 29:iMarital
	ReadCatValNum[29] =  2;		// 30:iMay75880
	ReadCatValNum[31] =  5;		// 32:iMillitary
	ReadCatValNum[32] =  3;		// 33:iMobility
	ReadCatValNum[33] =  3;		// 34:iMobillim
	ReadCatValNum[35] =  2;		// 36:iOthrserv
	ReadCatValNum[36] =  3;		// 37:iPerscare
	ReadCatValNum[47] =  2;		// 48:iRownchild
	ReadCatValNum[50] =  2;		// 51:iRrelchld
	ReadCatValNum[54] =  2;		// 55:iSept80
	ReadCatValNum[55] =  2;		// 56:iSex
	ReadCatValNum[60] =  2;		// 61:iVietnam
	ReadCatValNum[64] =  2;		// 65:iWWII
	ReadCatValNum[65] = 18;		// 66:iYearsch

	ReadUSCensus(ini.uscen_file, ini.read_cat, ini.read_cat_num, ini.catid_file);

	Exit(0);
}