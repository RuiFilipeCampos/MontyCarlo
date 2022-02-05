import React, { ReactNode } from 'react'
import Link from 'next/link'
import Head from 'next/head'
import { Text, Flex, VStack } from '@chakra-ui/layout'

type Props = {
  children: ReactNode
  title?: string
}

const Layout = ({ children, title = 'This is the default title' }: Props) => (
  <Flex minH="100vh" minW="100vw" bg="white" >
    <Head>
      <title>{title}</title>
      <meta charSet="utf-8" />
      <meta name="viewport" content="initial-scale=1.0, width=device-width" />
    </Head>
    {children}
  </Flex>
)

export default Layout
