////////////////////////////////////////////////////////////////////////////////
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
////////////////////////////////////////////////////////////////////////////////

#include "Factory.hpp"

#include <gtest/gtest.h>

namespace Aurora
{
    /// Helper class for FactoryTest. This is the base class for Product1 and
    /// Product2.
    class BaseProduct
    {
    };

    /// Helper class for FactoryTest. This class is derived from BaseProduct.
    class Product1 : public BaseProduct
    {
    };


    /// Helper class for FactoryTest. This class is derived from BaseProduct.
    class Product2 : public BaseProduct
    {
    };

    /// Test the factory template
    TEST(Factory, General)
    {
        Factory<BaseProduct> factory;

        // Register two products
        EXPECT_NO_THROW(factory.add<Product1>("Product1"));
        // Registering the same product twice should throw
        EXPECT_THROW(factory.add<Product2>("Product1"), EInvalidParameter);
        // Register Product2 now with its correct name
        EXPECT_NO_THROW(factory.add<Product2>("Product2"));

        // Pointer to an instance of the product
        Factory<BaseProduct>::BaseClassPtr p;
        // Create an instance of "Product1". As "Product1" has been correctly
        // registered. So this should not throw.
        EXPECT_NO_THROW(p = factory.create("Product1"));
        // The instance should not be zero
        EXPECT_TRUE(p != nullptr);
        // We have not registered "Product3", so this should throw
        EXPECT_THROW(p = factory.create("Product3"), EInvalidParameter);
    }
}
