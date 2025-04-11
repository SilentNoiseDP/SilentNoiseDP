#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <iostream>
#include <map>
using namespace std;

#define INI_PATH		".\\ReadFoursquare.ini"
#define CAT_NUM			415
#define USER_NUM	266909

struct iniParameter{
	char foursquare_dir[129];
	char catid_file[129];
	double v_min;
	double v_max;
	double h_min;
	double h_max;
};

struct location_type{
	double x;
	double y;
	char poi_name[129];
};

struct city_type{
	char name[129];
	double x;
	double y;
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
	char s[129];

	GetPrivateProfileString("PATH", "FOURSQUARE_FILE", "default", ret.foursquare_dir, sizeof(ret.foursquare_dir), INI_PATH);
	GetPrivateProfileString("PATH", "CATID_FILE", "default", ret.catid_file, sizeof(ret.catid_file), INI_PATH);

	GetPrivateProfileString("PARAMETER", "V_MIN", "default", s, sizeof(s), INI_PATH);
	ret.v_min = atof(s);
	GetPrivateProfileString("PARAMETER", "V_MAX", "default", s, sizeof(s), INI_PATH);
	ret.v_max = atof(s);
	GetPrivateProfileString("PARAMETER", "H_MIN", "default", s, sizeof(s), INI_PATH);
	ret.h_min = atof(s);
	GetPrivateProfileString("PARAMETER", "H_MAX", "default", s, sizeof(s), INI_PATH);
	ret.h_max = atof(s);

	return ret;
}

void ReadFoursquare(struct iniParameter ini){
	map<string, struct location_type> poi;
	string key;
	int userid;
	double user_x, user_y;
	char user_poi_name[129];
	char tim[1025];
	int user_cat;
	struct location_type loc;
	double min_dis, dis;
	char filename[1025];
	char s[1025];
	int i, j;
	char *tok, *ctx;
	FILE *fp, *fp2;

	sprintf_s(filename, "%s\\dataset_TIST2015_POIs_fix.txt", ini.foursquare_dir);
	FileOpen(&fp, filename, "r");

	i = 0;
	while (fgets(s, 1024, fp) != NULL){
		i++;

		tok = strtok_s(s, "\t", &ctx);
		key = string(tok);

		tok = strtok_s(NULL, "\t", &ctx);
		loc.y = atof(tok);

		tok = strtok_s(NULL, "\t", &ctx);
		loc.x = atof(tok);

		tok = strtok_s(NULL, "\t", &ctx);
		strcpy_s(loc.poi_name, tok);
		if (loc.poi_name[strlen(loc.poi_name) - 1] == '\n') loc.poi_name[strlen(loc.poi_name) - 1] = '\0';

		poi[key] = loc;
	}
	printf("POI File: %d lines\n", i);

	fclose(fp);

	sprintf_s(filename, "%s\\dataset_TIST2015_Checkins.txt", ini.foursquare_dir);
	FileOpen(&fp, filename, "r");

	FileOpen(&fp2, ini.catid_file, "w");

	i = 0;
	while (fgets(s, 1024, fp) != NULL){
		i++;

		tok = strtok_s(s, "\t", &ctx);
		userid = atoi(tok);

		tok = strtok_s(NULL, "\t", &ctx);
		key = string(tok);

		user_y = poi[key].y;

		user_x = poi[key].x;

		strcpy_s(user_poi_name, poi[key].poi_name);

		tok = strtok_s(NULL, "\t", &ctx);
		strcpy_s(tim, tok);

		if (userid <= 0 || userid > USER_NUM){
			printf("line:%d (user id:%d)\n", i, userid);
			Exit(-1);
		}

		if (user_y < ini.v_min || user_y > ini.v_max || user_x < ini.h_min || user_x > ini.h_max) continue;

		fprintf(fp2, "%d,%s,%f,%f,%s\n", userid, tim, user_x, user_y, user_poi_name);
	}

	fclose(fp);
	fclose(fp2);
}

void main(){
	struct iniParameter ini = readIniParameter();

	ReadFoursquare(ini);

	Exit(0);
}