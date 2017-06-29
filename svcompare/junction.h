/*
 *********************************************
 * getsv.h                                   *
 *                                           *
 *  Created on:   2011-5-24    v0.1.0        * 
 *  Modified on:  2013-1-10    v0.2.0        * 
 *  seeksv version:      0.8.0               * 
 *  Author:       qiukl                      * 
 *                                           * 
 *********************************************
 */

#pragma once

using namespace std;



class Junction
{
public:
	string up_chr;
	int up_pos;
	char up_strand;
	string down_chr;
	int down_pos;
	char down_strand;
	Junction() { }
	Junction(const string u_c, const int u_p, const char u_s, const string d_c, const int d_p, const char d_s) 
	{
		up_chr = u_c;
		up_pos = u_p;
		up_strand = u_s;
		down_chr = d_c;
		down_pos = d_p;
		down_strand = d_s;
	}
	void set_value(const string &u_c, const int &u_p, const char &u_s, const string &d_c, const int& d_p, const char& d_s) 
	{
		up_chr = u_c;
		up_pos = u_p;
		up_strand = u_s;
		down_chr = d_c;
		down_pos = d_p;
		down_strand = d_s;
	}
	bool operator< (const Junction& junction) const
	{
		if (up_chr < junction.up_chr)
			return 1;
		else if (up_chr > junction.up_chr)
			return 0;
		else
		{
			if (down_chr < junction.down_chr)
				return 1;
			else if (down_chr > junction.down_chr)
				return 0;
			else
			{
				if (up_strand < junction.up_strand)
					return 1;
				else if (up_strand > junction.up_strand)
					return 0;
				else
				{
					if (down_strand < junction.down_strand)
						return 1;
					else if (down_strand > junction.down_strand)
						return 0;
					else
					{
						// same chr , same strand , both upstream breakend and downstream breakend
						if (up_pos < junction.up_pos)
							return 1;
						else if (up_pos > junction.up_pos)
							return 0;
						else if (down_pos < junction.down_pos)
							return 1;
						else return 0;
					}
				}
			}
		}
	}
	friend ostream& operator<< (ostream &os, Junction &nj);
};

ostream& operator<< (ostream &os, Junction &nj)
{
	os << nj.up_chr << '\t' << nj.up_pos << '\t' << nj.up_strand << '\t'
       << nj.down_chr << '\t' << nj.down_pos << '\t' << nj.down_strand;
    return os;
}
