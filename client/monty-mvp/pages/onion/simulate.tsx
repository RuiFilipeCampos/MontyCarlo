import type { NextPage } from "next";
import React from "react";
import * as ch from "@chakra-ui/react";
import { NextRouter, useRouter } from "next/router";




const Onion: NextPage = (): JSX.Element => {
  // All states for handling the form information.
  return <>
  <ch.Flex w="100vw" h="100vh">
    <ch.HStack w="full" h="full" align="center" >
      <ch.VStack w="full" h="full" align="center" >
        <ch.Spacer />
        <ch.Spinner size='xl' />
        <ch.Text>Getting your Monty Carlo instance... ;)</ch.Text>
        <ch.Spacer />
      </ch.VStack>
    </ch.HStack>
  </ch.Flex>
</>
};

export default Onion;
