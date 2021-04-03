struct ShpHeader
{
	int file_code; //should be the integer 9994 in big endian byte order
	char blank[20]; //unused space
	unsigned length; //integer length in words (2 bytes) big endian byte order
	int version;
	int shape_type;
	double x_min;
	double y_min;
	double x_max;
	double y_max;
	double z_min;
	double z_max;
	double m_min;
	double m_max;
};

struct ShpRecHeader
{
	unsigned record_number;
	unsigned content_length;
};

struct Point
{
	double x_coordinate;
	double y_coordinate;
};

struct band
{
	double min;
	int numseg;
	int *seg;
};

struct ShpPointRecord
{
	int shape_type;
	double x;
	double y;
};

struct ShpPolygonRecord
{
	int shape_type;
	double x_min;
	double y_min;
	double x_max;
	double y_max;
	int num_parts;
	int num_points;
	int *parts; //index into points of the first point in each part.
	Point *points;
	int num_bands;
	band *ybands;
};

int read_file_header(FILE *in, ShpHeader *head);
int output_file_header(FILE *out, ShpHeader *head);
int read_record_header(FILE *in, ShpRecHeader *rechead);
int output_record_header(FILE *out, ShpRecHeader *rechead);
void read_points(char *filename, ShpHeader *head, ShpRecHeader **recheads, ShpPointRecord **records, int *numrec);
void cleanup_points(int numrec, ShpRecHeader *recheads, ShpPointRecord *records);
int read_poly_record(FILE *in, ShpPolygonRecord *rec);
int output_poly_record(FILE *out, ShpPolygonRecord *rec);
void read_polygons(char *filename, ShpHeader *head, ShpRecHeader **recheads, ShpPolygonRecord **records, int *numrec);
void cleanup_polygons(int numrec, ShpRecHeader *recheads, ShpPolygonRecord *records);
int in_interior(Point *p, ShpPolygonRecord *poly);
int get_id(Point *p, ShpPolygonRecord *records, int numrec);
