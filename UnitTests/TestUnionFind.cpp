#include "CppUnitTest.h"

#include "SignedUnionFind.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTests
{		
	TEST_CLASS(TestUnionFind)
	{
	public:
		template<typename TUnionFind>
		void testInitBody(TUnionFind& uf)
		{
			for (int i = 0; i < 5; ++i)
				uf.addItem(false);

			for (SignedUnionFind<true>::index_t i = 0; i < 5; ++i)
				Assert::AreEqual(i, uf.getRepresentative(i), L"Parent pointer of unmerged leaf does not point to itself.");
		}

		TEST_METHOD(TestInitWithSign)
		{
			SignedUnionFind<true> uf;

			testInitBody(uf);
		}

		TEST_METHOD(TestInitWithoutSign)
		{
			SignedUnionFind<false> uf;

			testInitBody(uf);
		}

		template<typename TUnionFind>
		void testSimpleMergeBody(TUnionFind& uf)
		{			
			uf.addItem(false);
			uf.addItem(false);

			uf.merge(0, 1);

			Assert::AreEqual(uf.getRepresentative(0),
				uf.getRepresentative(1), L"Representatives are not equal after merge.");
		}

		
		TEST_METHOD(TestSimpleMergeWithoutSign)
		{
			SignedUnionFind<false> uf;
			
			testSimpleMergeBody(uf);
		}

		TEST_METHOD(TestSimpleMergeWithSign)
		{
			SignedUnionFind<true> uf;

			testSimpleMergeBody(uf);
		}

		template <typename TUnionFind>
		void testMultipleMergeBody(TUnionFind& uf)
		{
			for (int i = 0; i < 8; ++i)
				uf.addItem(false);

			uf.merge(0, 1);
			uf.merge(2, 3);
			uf.merge(4, 5);
			uf.merge(6, 7);

			Assert::AreEqual(uf.getRepresentative(0), uf.getRepresentative(1));
			Assert::AreEqual(uf.getRepresentative(2), uf.getRepresentative(3));
			Assert::AreEqual(uf.getRepresentative(4), uf.getRepresentative(5));
			Assert::AreEqual(uf.getRepresentative(6), uf.getRepresentative(7));

			Assert::AreNotEqual(uf.getRepresentative(0), uf.getRepresentative(2));

			uf.merge(0, 2);
			uf.merge(5, 7);

			Assert::AreEqual(uf.getRepresentative(1), uf.getRepresentative(3));
			Assert::AreEqual(uf.getRepresentative(4), uf.getRepresentative(6));

			uf.merge(1, 6);
			for (int i = 0; i < 8; ++i)
				for (int j = 0; j < 8; ++j)
					Assert::AreEqual(uf.getRepresentative(i), uf.getRepresentative(j));
		}

		TEST_METHOD(TestMultipleMergeWithoutSign)
		{
			SignedUnionFind<false> uf;
			
			testMultipleMergeBody(uf);
		}

		TEST_METHOD(TestMultipleMergeWithSign)
		{
			SignedUnionFind<true> uf;

			testMultipleMergeBody(uf);
		}

		TEST_METHOD(TestSimpleSign)
		{
			SignedUnionFind<true> uf;
			for (int i = 0; i < 4; ++i)
				uf.addItem(i % 2 == 0);			

			for (int i = 0; i < 4; ++i)
				Assert::AreEqual(i % 2 == 0, uf.getSign(i));
		}

		TEST_METHOD(TestInflatedSign)
		{
			SignedUnionFind<true> uf;
			for (int i = 0; i < 4; ++i)
				uf.addItem(i % 2 == 0);

			uf.merge(0, 1);
			uf.merge(3, 2);
			uf.merge(1, 3);

			for (int i = 0; i < 4; ++i)
				Assert::AreEqual(i % 2 == 0, uf.getSign(i));
		}

		TEST_METHOD(TestFlipMergedSign)
		{
			SignedUnionFind<true> uf;
			for (int i = 0; i < 4; ++i)
				uf.addItem(i % 2 == 0);
			uf.addItem(false);

			uf.merge(0, 1);
			uf.merge(3, 2);
			uf.merge(1, 3);

			for (int j = 0; j < 4; ++j)
			{
				uf.flipSign(j); //changing any sign should affect all other nodes except 4
				for (int i = 0; i < 4; ++i)
					Assert::AreEqual(i % 2 != j % 2, uf.getSign(i));
				Assert::AreEqual(false, uf.getSign(4));
			}
		}
	};
}