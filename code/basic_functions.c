char* concat(const char *s1, const char *s2){
	char *result = (char*)malloc(strlen(s1)+strlen(s2)+1);
	strcpy(result, s1);
	strcat(result, s2);
	return result;
}

char* concat(const char *s1, const char *s2, const char *s3){
	char *result = (char*)malloc(strlen(s1)+strlen(s2)+strlen(s3)+1);
	strcpy(result, s1);
	strcat(result, s2);
	strcat(result, s3);
	return result;
}

char* concat(const char *s1, const char *s2, const char *s3, const char *s4){
	char *result = (char*)malloc(strlen(s1)+strlen(s2)+strlen(s3)+strlen(s4)+1);
	strcpy(result, s1);
	strcat(result, s2);
	strcat(result, s3);
	strcat(result, s4);
	return result;
}

char* concat(const char *s1, const char *s2, const char *s3, const char *s4, const char *s5){
	char *result = (char*)malloc(strlen(s1)+strlen(s2)+strlen(s3)+strlen(s4)+strlen(s5)+1);
	strcpy(result, s1);
	strcat(result, s2);
	strcat(result, s3);
	strcat(result, s4);
	strcat(result, s5);
	return result;
}
