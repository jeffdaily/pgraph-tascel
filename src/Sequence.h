#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <string>

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
        explicit Sequence();

        /** shallow copy; new instance does not own the data */
        explicit Sequence(const Sequence &seq);

        /** from a std::string */
        explicit Sequence(const string &seq);

        /** from character sequence; new instance does not own the data */
        Sequence(const char *seq, size_t n);

        /** delete data, if owner */
        virtual ~Sequence();

    private:
        bool is_owner;      /**< whether this instance should delete the data */
        const char *data;   /**< character sequence */
        size_t size;        /**< length of character sequence */
};

#endif /* SEQUENCE_H_ */
