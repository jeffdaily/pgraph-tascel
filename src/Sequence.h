#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <iostream>
#include <string>

using std::ostream;
using std::string;


/**
 * A biological sequence.
 *
 * Similar to std::string, but immutable and the memory potentially may not be
 * owned by any instance but rather as part of a larger SequenceDatabase
 * instance (or some other data source).
 */
class Sequence
{
    public:
        /** default */
        Sequence();

        /** shallow copy; new instance does not own the data */
        Sequence(const Sequence &seq);

        /** from a std::string */
        Sequence(const string &seq);

        /** from character sequence; new instance does not own the data */
        Sequence(const char *seq);

        /** from character sequence; new instance does not own the data */
        Sequence(const char *seq, size_t len);

        /** delete data, if owner */
        virtual ~Sequence();

        /** affine gap align, returning score, #matches, and align length */
        void align(const Sequence &that, int &score, int &ndig, int &alen);

        /** overload of string cast */
        operator string() const;

        friend ostream& operator << (ostream &os, const Sequence &s);

    private:
        bool is_owner;      /**< whether this instance should delete the data */
        const char *data;   /**< character sequence */
        size_t size;        /**< length of character sequence */
};


#endif /* SEQUENCE_H_ */
