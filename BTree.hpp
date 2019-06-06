//
// Created by 郑文鑫 on 2019-03-09.
//

#include "utility.hpp"
#include <functional>
#include <cstddef>
#include "exception.hpp"
#include <cstdio>
namespace sjtu {
    template <class Key, class Value, class Compare = std::less<Key> >
    class BTree {
	pubilc:
		typedef pair<const Key, Value> value_type;
		typedef size_t off_t;
		typedef pair<Key, off_t> Keynchil;
	private:
	
		constexpr off_t node_size = 4096;
		constexpr off_t key_size = sizeof(Key);
		constexpr off_t value_size = sizeof(Value);
		constexpr off_t L = (node_size -sizeof(int)-sizeof(size_t)*4) / sizeof(key_size + value_size) - 1;
		constexpr off_t M = (node_size - sizeof(bool) - sizeof(int) - sizeof(size_t) * 2) / (sizeof(Keynchil)) - 1;

		struct info_node {
			off_t root=0;
			off_t head=0;
			off_t tail=0;
			off_t total_size=0;
			off_t eof=0;

		};
		
		struct par_node {
			bool type;//如果其child是leaf则为1，否则为0
			off_t offset;
			off_t par;
			int num;
			Keynchil keynchil[M+1];
		};
		struct leaf_node {
			off_t offset;
			off_t par;
			int num;
			off_t prev, next;
			value_type data[L + 1];
		};

		FILE* fp=NULL;
		bool file_open=0;
		info_node info;


		void read_node(void *buffer, off_t offset,size_t buff_size)
		{
			fseek(fp, offset, 0);
			fread(buffer,buff_size , 1, fp);
		}


		void write_node(void* buffer, off_t offset,size_t buff_size)
		{
			fseek(fp, offset, 0);
			fwrite(buffer, buff_size, 1, fp);
		}

		void build_tree()
		{
			info.total_size = 0;
			
			par_node p;
			leaf_node leaf;
			p.type = 1;
			p.offset = sizeof(info_node);
			p.par = 0;
			p.num=0;
			leaf.offset = sizeof(info_node) + sizeof(par_node);
			info.eof = leaf.offset + sizeof(leaf_node);
			info.root = p.offset;
			info.head = info.tail = leaf.offset;
			p.keynchil[0].second = leaf.offset;
			leaf.par = p.offset;
			leaf.num = 0;
			p.prev = p.next = leaf.prev = leaf.next = 0;
			write_node(&info, 0, sizeof(info_node));
			write_node(&p, p.offset, sizeof(par_node));
			write_node(&leaf, leaf.offset, sizeof(leaf_node));
		}

		void add_to_par(off_t chil, Key key, par_node &par)
		{
			int i;
			for (i = 0; i < par.num; ++i)
				if (key < par.keynchil[i].first)
					break;
			for (int j = par.num; j > i; --j)
				par.keynchil[j] = par.keynchil[j - 1];
			++par.num;
			par.keynchil[i].first = key;
			par.keynchil[i].second = chil;
			if (par.num > M)
				split_par(par);
			else write_node(par, par.offset, sizeof(par_node));
		}

		void split_leaf(leaf_node &oldleaf)
		{
			leaf_node newleaf;
			newleaf.num = oldleaf.num >> 1;
			oldleaf.num -= newleaf.num;
			newleaf.offset = info.eof;
			info.eof += sizeof(leaf_node);
			newleaf.par = oldlaef.par;
			newleaf.next = oldleaf.next;
			newleaf.prev = oldleaf.offset;
			oldleaf.next = newleaf.offset;
			if (newleaf.next)
			{
				leaf_node tmpleaf;
				read_node(&tmpleaf, newleaf.next, sizeof(leaf_node));
				tmpleaf.prev = newleaf.offset;
				write_node(&tmpleaf, tmpleaf.offset, sizeof(leaf_node));
			}
			else {
				info.tail = newleaf.offset;
			}
			for (int i = 0; i < newleaf.num; ++i)
			{
				newleaf.data[i] = oldleaf.data[i + oldleaf.num];
			}
			par_node p;
			read_node(&p, oldleaf.par, sizeof(par_node));
			write_node(&info, 0, sizeof(info_node));
			write_node(&oldleaf, oldleaf.offset, sizeof(leaf_node));
			write_node(&newleaf, newleaf.offset, sizeof(leaf_node));
			add_to_par(newleaf.offset, newleaf.data[0].first, p);
		}


		void split_par(par_node &p) {
			par_node newp;
			newp.num = p.num >> 1;
			p.num -= newp.num;
			newp.par = p.par;
			newp.type = p.type;
			newp.offset = info.eof;
			info.eof += sizeof(par_node);
			for (int i = 0; i < newp.num; ++i)
				newp.keynchil[i] = p.keynchil[i + p.num];
			if (newp.type)
			{
				leaf_node tmp;
				for (int i = 0; i < newp.num; ++i)
				{
					read_node(&tmp, newp.keynchil[i].second, sizeof(leaf_node));
					tmp.par = newp.offset;
					write_node(&tmp, tmp.offset, sizeof(leaf_node));
				}
			}
			else {
				par_node tmp;
				for (int i = 0; i < newp.num; ++i)
				{
					read_node(&tmp, newp.keynchil[i].second, sizeof(par_node));
					tmp.par = newp.offset;
					write_node(&tmp, tmp.offset, sizeof(par_node));
				}
			}
			if (info.root == p.offset)
			{
				par_node nroot;
				nroot.par = 0;
				nroot.num = 2;
				nroot.type = 0;
				nroot.offset = info.eof;
				info.eof += sizeof(par_node);
				info.root = nroot.offset;
				nroot.keynchil[0].first = p.keynchil[0].first;
				nroot.keynchil[1].first = newp.keynchil[0].first;
				nroot.keynchil[0].second = p.offset;
				nroot.keynchil[1].second = newp.offset;
				p.par = nroot.offset;
				newp.par = nroot.offset;
				write_node(&p, p.offset, sizeof(par_node));
				write_node(&newp, newp.offset, sizeof(par_node));
				write_node(&nroot, nroot.offset, sizeof(par_node));
				write_node(&info, 0, sizeof(par_node));
			}
			else {
				write_node(&p, p.offset, sizeof(par_node));
				write_node(&newp, newp.offset, sizeof(par_node));
				write_node(&info, 0, sizeof(par_node));
				par_node par;
				read_node(&par, p.par, sizeof(par_node));
				add_to_par(newp.offset, newp.keynchil[0].first, par);
			}
		}

		bool leaf_borrow_next(leaf_node leaf)
		{
			if (!leaf.next) return 0;
			leaf_node tmp;
			read_node(&tmp, leaf.next, sizeof(leaf_node));
			if (leaf.par != tmp.par || tmp.num <= L / 2)
				return 0;
			leaf.data[leaf.num++] = tmp.data[0];
			--tmp.num;
			for (int i = 0; i < tmp.num; ++i)
				tmp.data[i] = tmp.data[i + 1];
			par_node par;
			read_node(&par, leaf.par, sizeof(par_node));
			for (int i=0;i<par.num;++i)
				if (par.keynchil[i].first == leaf.data[leaf.num - 1].first)
				{
					par.keynchil[i].first == tmp.data[0].first;
					break;
				}
			write_node(&par, par.offset, sizeof(par_node));
			write_node(&leaf, leaf.offset, sizeof(leaf_node));
			write_node(&tmp, tmp.offset, sizeof(leaf_node));
			return 1;
		}
			
		bool leaf_borrow_prev(leaf_node leaf)
		{
			if (!leaf.prev) return 0;
			leaf_node tmp;
			read_node(&tmp, leaf.prev, sizeof(leaf_node));
			if (leaf.par != tmp.par || tmp.num <= L / 2)
				return 0;
			for (int i = leaf.num; i >0;--i)
				leaf.data[i] = leaf.data[i - 1];
			leaf.data[0] = tmp.data[tmp.num - 1];
			++leaf.num;
			--tmp.num;

			par_node par;
			read_node(&par, leaf.par, sizeof(par_node));
			for (int i = 0; i < par.num; ++i)
				if (par.keynchil[i].first == leaf.data[1].first)
				{
					par.keynchil[i].first == leaf.data[0].first;
					break;
				}
			write_node(&par, par.offset, sizeof(par_node));
			write_node(&leaf, leaf.offset, sizeof(leaf_node));
			write_node(&tmp, tmp.offset, sizeof(leaf_node));
			return 1;
		}

		bool merge_leaf_next(leaf_node leaf)
		{
			if (!leaf.next) return 0;
			leaf_node tmp;
			read_node(&tmp, leaf.next, sizeof(leaf_node));
			if (leaf.par != tmp.par)
				return 0;
			leaf.next = tmp.next;
			if (info.tail == tmp.offset)
				info.tail == leaf.offset;
			else {
				leaf_node nex;
				read_node(&nex, leaf.next; sizeof(leaf_node));
				nex.prev = leaf.offset;
				write_node(&nex, nex.offset, sizeof(leaf_node));
			}
			for (int i = 0; i < tmp.num; ++i)
				leaf.data[i + leaf.num] = tmp.data[i];
			leaf.num += tmp.num;
			par_node par;
			read_node(par, leaf.par, sizeof(par_node));
			int i;
			for (i = 0; i < par.num; ++i)
				if (par.keynchil[i].first == tmp.data[0].first)
					break;
			for (int j = i; j < par.num - 1; ++j)
				par.keynchil[j] = par.keynchil[j + 1];
			--par.num;
			write_node(&par, par.offset, sizeof(par_node));
			write_node(&leaf, leaf.offset, sizeof(leaf_node));
			return 1;
		}

		bool merge_leaf_prev(leaf_node leaf)
		{
			if (!leaf.prev) return 0;
			leaf_node tmp;
			read_node(&tmp, leaf.prev, sizeof(leaf_node));
			if (leaf.par != tmp.par)
				return 0;
			tmp.next = leaf.next;
			if (info.tail == leaf.offset)
				info.tail == tmp.offset;
			else {
				leaf_node nex;
				read_node(&nex, leaf.next; sizeof(leaf_node));
				nex.prev = tmp.offset;
				write_node(&nex, nex.offset, sizeof(leaf_node));
			}
			for (int i = 0; i < leaf.num; ++i)
				tmp.data[i + tmp.num] = leaf.data[i];
			tmp.num += leaf.num;
			par_node par;
			read_node(par, leaf.par, sizeof(par_node));
			int i;
			for (i = 0; i < par.num; ++i)
				if (par.keynchil[i].first == leaf.data[0].first)
					break;
			for (int j = i; j < par.num - 1; ++j)
				par.keynchil[j] = par.keynchil[j + 1];
			--par.num;
			write_node(&par, par.offset, sizeof(par_node));
			write_node(&tmp, tmp.offset, sizeof(leaf_node));
			return 1;


		}

		bool par_borrow_next(par_node p)
		{
			if (!p.par) return 0;
			par_node par;
			read_node(&par, p.par, sizeof(par_node));
			int i;
			for (i = 0; i < par.num; ++i)
				if (par.keynchil[i].first == p.keynchil[0].first)
				    break;
			if (i + 1 >= par.num) return 0;
			par_node tmp;
			read_node(&tmp, par.keynchil[i+1].second, sizeof(par_node));
			if (tmp.num <= M / 2)
				return 0;
			p.keynchil[p.num++] = tmp.keynchil[0];
			--tmp.num;
			for (int i = 0; i < tmp.num; ++i)
				tmp.keynchil[i] = tmp.keynchil[i + 1];
			par.keynchil[i + 1].first = tmp.keynchil[0].first;
			if (p.type == 1)
			{
				leaf_node ch;
				read_node(&ch, p.keynchil[p.num - 1].second, sizeof(leaf_node));
				ch.par = p.offset;
				write_node(&ch, ch.offset, sizeof(leaf_node));
			}
			else {

				par_node ch;
				read_node(&ch, p.keynchil[p.num - 1].second, sizeof(par_node));
				ch.par = p.offset;
				write_node(&ch, ch.offset, sizeof(par_node));
			}
			
			write_node(&par, par.offset, sizeof(par_node));
			write_node(&p, p.offset, sizeof(par_node));
			write_node(&tmp, tmp.offset, sizeof(par_node));
			return 1;


		}

		bool par_borrow_prev(par_node p)
		{
			if (!p.par) return 0;
			par_node par;
			read_node(&par, p.par, sizeof(par_node));
			int i;
			for (i = 0; i < par.num; ++i)
				if (par.keynchil[i].first == p.keynchil[0].first)
					break;
			if (i - 1 < 0) return 0;
			par_node tmp;
			read_node(&tmp, par.keynchil[i - 1].second, sizeof(par_node));
			if (tmp.num <= M / 2)
				return 0;
			for (int j = p.num; j > 0; ++j)
				p.keynchil[j] = p.keynchil[j - 1];
			p.keynchil[0] = tmp.keynchil[tmp.num-1];
			++p.num;
			--tmp.num;
			par.keynchil[i].first = p.keynchil[0].first;
			if (p.type == 1)
			{
				leaf_node ch;
				read_node(&ch, p.keynchil[0].second, sizeof(leaf_node));
				ch.par = p.offset;
				write_node(&ch, ch.offset, sizeof(leaf_node));
			}
			else {

				par_node ch;
				read_node(&ch, p.keynchil[0].second, sizeof(par_node));
				ch.par = p.offset;
				write_node(&ch, ch.offset, sizeof(par_node));
			}

			write_node(&par, par.offset, sizeof(par_node));
			write_node(&p, p.offset, sizeof(par_node));
			write_node(&tmp, tmp.offset, sizeof(par_node));
			return 1;
		}

		bool merge_par_next(par_node p)
		{
			if (!p.par) return 0;
			par_node par;
			read_node(&par, p.par, sizeof(par_node));
			int i;
			for (i = 0; i < par.num; ++i)
				if (par.keynchil[i].first == p.keynchil[0].first)
					break;
			if (i >= par.num - 1)
				return 0;

			par_node tmp;
			read_node(&tmp, par.keynchil[i+1].second, sizeof(par_node));
			
			for (int i = 0; i < tmp.num; ++i)
				p.data[i + leaf.num] = tmp.data[i];
			p.num += tmp.num;
			for (int j = i + 1; j < par.num - 1; ++j)
				par.keynchil[j] = par.keynchil[j + 1];
			--par.num;
			write_node(&par, par.offset, sizeof(par_node));
			write_node(&leaf, leaf.offset, sizeof(leaf_node));
			return 1;
		}

		bool merge_par_prev(par_node p)
		{
			if (!p.par) return 0;
			par_node par;
			read_node(&par, p.par, sizeof(par_node));
			int i;
			for (i = 0; i < par.num; ++i)
				if (par.keynchil[i].first == p.keynchil[0].first)
					break;
			if (i == 0)
				return 0;

			par_node tmp;
			read_node(&tmp, par.keynchil[i - 1].second, sizeof(par_node));

			for (int i = 0; i < p.num; ++i)
				tmp.data[i + tmp.num] = p.data[i];
			tmp.num += p.num;
			for (int j = i; j < par.num - 1; ++j)
				par.keynchil[j] = par.keynchil[j + 1];
			--par.num;
			write_node(&par, par.offset, sizeof(par_node));
			write_node(&tmp, tmp.offset, sizeof(leaf_node));
			return 1;
		}

		off_t find_pos(Key key, off_t paroff)
		{
			par_node par;
			read_node(&par, paroff, sizeof(par_node));



		}

    public:
        
        
        
        class const_iterator;
        class iterator {
        private:
            // Your private members go here
        public:
            bool modify(const Value& value){
                
            }
            iterator() {
                // TODO Default Constructor
            }
            iterator(const iterator& other) {
                // TODO Copy Constructor
            }
            // Return a new iterator which points to the n-next elements
            iterator operator++(int) {
                // Todo iterator++
            }
            iterator& operator++() {
                // Todo ++iterator
            }
            iterator operator--(int) {
                // Todo iterator--
            }
            iterator& operator--() {
                // Todo --iterator
            }
            // Overloaded of operator '==' and '!='
            // Check whether the iterators are same
            bool operator==(const iterator& rhs) const {
                // Todo operator ==
            }
            bool operator==(const const_iterator& rhs) const {
                // Todo operator ==
            }
            bool operator!=(const iterator& rhs) const {
                // Todo operator !=
            }
            bool operator!=(const const_iterator& rhs) const {
                // Todo operator !=
            }
        };
        class const_iterator {
            // it should has similar member method as iterator.
            //  and it should be able to construct from an iterator.
        private:
            // Your private members go here
        public:
            const_iterator() {
                // TODO
            }
            const_iterator(const const_iterator& other) {
                // TODO
            }
            const_iterator(const iterator& other) {
                // TODO
            }
            // And other methods in iterator, please fill by yourself.
        };
        // Default Constructor and Copy Constructor
        BTree() {
			if (file_open) return;
			fp = fopen("bpt.dat", "rb+");
			if (!fp)
			{
				fp = fopen("bpt.dat", "wb+");
				build_tree();
			}
			else {
				fseek(fp, 0, 0);
				fread(info, 1, sizeof(info_node), 1, fp);
			}
			file_open = 1;
        }
        BTree(const BTree& other) {
			fp = fopen("bpt.dat", "rb+");
			info.eof = other.info.eof;
			info.root = other.info.root;
			info.head = other.info.head;
			info.tail = other.info.tail;
			info.total_size = other.total_size;
			file_open = 1;
        }
        BTree& operator=(const BTree& other) {
			fp = fopen("bpt.dat", "rb+");
			info.eof = other.info.eof;
			info.root = other.info.root;
			info.head = other.info.head;
			info.tail = other.info.tail;
			info.total_size = other.total_size;
			file_open = 1;
			return *this;
        }
        ~BTree() {
			if (file_open)
			{
				file_open = 0;
				fclose(fp);
			}
        }
        // Insert: Insert certain Key-Value into the database
        // Return a pair, the first of the pair is the iterator point to the new
        // element, the second of the pair is Success if it is successfully inserted
        pair<iterator, OperationResult> insert(const Key& key, const Value& value) {
			off_t leafoff = find_pos(key, info.root);


        }
        // Erase: Erase the Key-Value
        // Return Success if it is successfully erased
        // Return Fail if the key doesn't exist in the database
        OperationResult erase(const Key& key) {
            // TODO erase function
            return Fail;  // If you can't finish erase part, just remaining here.
        }
        // Return a iterator to the beginning
        iterator begin() {}
        const_iterator cbegin() const {}
        // Return a iterator to the end(the next element after the last)
        iterator end() {}
        const_iterator cend() const {}
        // Check whether this BTree is empty
        bool empty() const {}
        // Return the number of <K,V> pairs
        size_t size() const {}
        // Clear the BTree
        void clear() {}
        // Return the value refer to the Key(key)
        Value at(const Key& key){
        }
        /**
         * Returns the number of elements with key
         *   that compares equivalent to the specified argument,
         * The default method of check the equivalence is !(a < b || b > a)
         */
        size_t count(const Key& key) const {}
        /**
         * Finds an element with key equivalent to key.
         * key value of the element to search for.
         * Iterator to an element with key equivalent to key.
         *   If no such element is found, past-the-end (see end()) iterator is
         * returned.
         */
        iterator find(const Key& key) {}
        const_iterator find(const Key& key) const {}
    };
}  // namespace sjtu

